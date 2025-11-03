#!/usr/bin/env python3

"""
Production-Ready MD Trajectory Processor

This script extracts molecular geometries and velocities from large-scale
MD (Molecular Dynamics) trajectory files.

It is designed for performance, reliability, and user-friendliness, incorporating:
- Parallel (multiprocessing) frame counting
- Threaded, simultaneous data extraction for I/O efficiency
- TQDM progress bars
- Detailed logging and summary reporting
- Multiprocessing-based data validation
- Dynamic resource management and reporting
- Graceful interrupt handling

Requirements:
- tqdm
- psutil
"""

import os
import sys
import argparse
import multiprocessing
import logging
import time
import csv
import random
import psutil
from tqdm import tqdm
from threading import Thread, Event
from datetime import datetime

# --- Configuration ---
# Use a specific marker to find frames in the coordinate file
FRAME_MARKER = "FINAL HEAT OF FORMATION ="
DEFAULT_CHUNK_MB = 100  # Default chunk size for reading
LARGE_FILE_FRAMES = 50000  # Threshold to trigger dynamic resource scaling
VALIDATION_SAMPLE_SIZE = 100  # Default number of frames to validate

# --- Banner ---
def print_banner(log_level=logging.INFO):
    """Prints a banner to the log."""
    banner = """
=================================================
     MD Trajectory Production Processor v2.0
=================================================
    """
    logging.log(log_level, banner)


# --- 1. Robust Frame Counting (Multiprocessing) ---

def _count_worker(args):
    """Internal worker for count_pattern_in_file."""
    file_path, pattern, start_byte, end_byte = args
    pattern_len = len(pattern)
    overlap = pattern_len - 1
    read_start = max(0, start_byte - overlap)
    read_size = (end_byte - read_start) + overlap
    count = 0
    
    try:
        with open(file_path, "r", errors="ignore") as f:
            f.seek(read_start)
            data = f.read(read_size)
            
            idx = data.find(pattern)
            while idx != -1:
                abs_pos = read_start + idx
                # Count if the pattern *starts* within our assigned chunk
                if start_byte <= abs_pos < end_byte:
                    count += 1
                idx = data.find(pattern, idx + 1)
        return count
    except Exception as e:
        logging.error(f"[CountWorker] Error: {e}")
        return 0

def count_pattern_in_file(file_path, pattern, num_workers):
    """
    Counts occurrences of a pattern in a large text file using multiple processes.
    
    Args:
        file_path (str): Path to the file.
        pattern (str): String pattern to search for.
        num_workers (int): Number of CPU cores to use.

    Returns:
        int: Total count of the pattern.
    """
    try:
        file_size = os.path.getsize(file_path)
    except FileNotFoundError:
        logging.error(f"File not found for counting: {file_path}")
        raise
    except Exception as e:
        logging.error(f"Could not get size of file {file_path}: {e}")
        raise

    if file_size == 0:
        logging.warning(f"File {file_path} is empty. Frame count is 0.")
        return 0

    # Calculate chunk ranges for each worker
    chunk_size = file_size // num_workers
    ranges = []
    for i in range(num_workers):
        start_byte = i * chunk_size
        end_byte = (i + 1) * chunk_size if i < num_workers - 1 else file_size
        ranges.append((file_path, pattern, start_byte, end_byte))

    logging.info(f"Starting frame count with {num_workers} workers...")
    
    total_count = 0
    with tqdm(total=num_workers, desc="Counting Frames", unit="worker", leave=False) as pbar:
        with multiprocessing.Pool(processes=num_workers) as pool:
            for count in pool.imap_unordered(_count_worker, ranges):
                total_count += count
                pbar.update(1)
                
    logging.info(f"Frame counting complete. Found {total_count} occurrences.")
    return total_count


# --- 2. Core Data Streaming Functions ---

def process_coord_chunk(chunk, leftover=""):
    """Processes a chunk of a coordinate .xyz file."""
    results = []
    lines = (leftover + chunk).splitlines()
    i = 0
    leftover_idx = len(lines)

    while i < len(lines):
        try:
            # Find the number of atoms
            num_atoms = int(lines[i].strip())
            frame_size = num_atoms + 2  # N_atoms + title_line + data_lines
        except (ValueError, IndexError):
            # Not a valid frame start, skip
            i += 1
            continue

        if i + frame_size > len(lines):  # Incomplete frame
            leftover_idx = i
            break
        
        frame_data = lines[i : i + frame_size]
        frame_results = []
        # Line 0 is N_atoms, Line 1 is title/comment
        for atom_line in frame_data[2:]:
            parts = atom_line.split()
            if len(parts) == 4:
                try:
                    atom, x, y, z = parts
                    frame_results.append((atom, float(x), float(y), float(z)))
                except ValueError:
                    logging.debug(f"Skipping malformed coord line: {atom_line}")
            else:
                 logging.debug(f"Skipping malformed coord line (wrong part count): {atom_line}")

        # Check if frame is complete before appending
        if len(frame_results) == num_atoms:
            results.append((num_atoms, lines[i+1], frame_results)) # (N, title, data)
        else:
            logging.warning(f"Frame starting at line index {i} reports {num_atoms} atoms but found {len(frame_results)}. Skipping frame.")
            
        i += frame_size
        leftover_idx = i

    final_leftover = "\n".join(lines[leftover_idx:])
    return results, final_leftover

def chunk_reader_coords(file_path, chunk_size, stop_event):
    """Generator that streams coordinate frames from a file."""
    frame_index = 0
    leftover = ""
    try:
        with open(file_path, "r") as f:
            while not stop_event.is_set():
                chunk = f.read(chunk_size)
                if not chunk:
                    if leftover: # Process final leftover data
                        results, _ = process_coord_chunk(chunk, leftover)
                        for frame in results:
                            frame_index += 1
                            yield frame_index, frame
                    break
                
                results, leftover = process_coord_chunk(chunk, leftover)
                for frame in results:
                    frame_index += 1
                    yield frame_index, frame
    except FileNotFoundError:
        logging.error(f"Coords file not found: {file_path}")
        stop_event.set() # Signal other threads to stop
    except Exception as e:
        logging.error(f"Error reading coords file {file_path}: {e}")
        stop_event.set()

def process_velocity_chunk(chunk, leftover="", atoms_per_frame=0):
    """Processes a chunk of a velocity .vel file."""
    results = []
    lines = (leftover + chunk).splitlines()
    frame = []
    in_velocity_block = False
    
    # Check if leftover had an open velocity block
    if "START_VELOC_BLOCK" in leftover:
        in_velocity_block = True
    
    # Use line numbers for better leftover handling
    leftover_idx = len(lines)

    for i, line in enumerate(lines):
        if "$VELOC" in line:
            in_velocity_block = True
            frame = [] # Start a new frame
            continue
        
        if "$ENDVELOC" in line:
            in_velocity_block = False
            if atoms_per_frame and len(frame) == atoms_per_frame:
                results.append(frame)
            elif not atoms_per_frame and frame:
                 results.append(frame) # Append if dynamic
            frame = []
            continue

        if not in_velocity_block:
            continue

        parts = line.split()
        if len(parts) == 3:
            try:
                frame.append(tuple(map(float, parts)))
            except ValueError:
                logging.debug(f"Skipping malformed velocity line: {line}")

        # If we know the atoms_per_frame, we can save completed frames
        # *before* $ENDVELOC, but this is less common for MOPAC.
        # Sticking to $ENDVELOC as the frame delimiter is safer.
        
        # We assume the *next* line might be part of the current frame,
        # so we only set leftover_idx when *outside* a block.
        leftover_idx = i + 1

    # Handle leftover
    if in_velocity_block:
        # If we are still in a block, the entire block is leftover
        # Find start of block in this chunk
        block_start_idx = -1
        for i, line in enumerate(lines):
             if "$VELOC" in line:
                 block_start_idx = i
        
        if block_start_idx != -1:
             final_leftover = "START_VELOC_BLOCK\n" + "\n".join(lines[block_start_idx:])
        else:
             # Block started in a *previous* chunk
             final_leftover = "START_VELOC_BLOCK\n" + "\n".join(lines)
    else:
        # Not in a block, no leftover
        final_leftover = ""

    return results, final_leftover

def chunk_reader_velocity(file_path, atoms_per_frame, chunk_size, stop_event):
    """Generator that streams velocity frames from a file."""
    frame_index = 0
    leftover = ""
    try:
        with open(file_path, "r") as f:
            while not stop_event.is_set():
                chunk = f.read(chunk_size)
                if not chunk:
                    if leftover:
                         results, _ = process_velocity_chunk(chunk, leftover, atoms_per_frame)
                         for frame in results:
                            frame_index += 1
                            yield frame_index, frame
                    break
                
                results, leftover = process_velocity_chunk(chunk, leftover, atoms_per_frame)
                for frame in results:
                    frame_index += 1
                    yield frame_index, frame
    except FileNotFoundError:
        logging.error(f"Velocity file not found: {file_path}")
        stop_event.set() # Signal other threads to stop
    except Exception as e:
        logging.error(f"Error reading velocity file {file_path}: {e}")
        stop_event.set()

def infer_atoms_from_velocity_file(vel_path, max_lines=50000):
    """
    Infers the number of atoms per velocity frame by scanning the first $VELOC block.
    """
    logging.info(f"Attempting to infer atoms_per_frame from {vel_path}...")
    in_block = False
    count = 0
    try:
        with open(vel_path, "r") as f:
            for i, line in enumerate(f):
                if i > max_lines:
                    logging.warning(f"Could not find velocity block in first {max_lines} lines.")
                    return None
                if "$VELOC" in line:
                    in_block = True
                    count = 0
                    continue
                if "$ENDVELOC" in line and in_block:
                    logging.info(f"Inferred {count} atoms per frame.")
                    return count
                if in_block:
                    parts = line.split()
                    if len(parts) == 3:
                        count += 1
    except Exception as e:
        logging.error(f"Error during velocity atom inference: {e}")
        return None
    return None

def infer_atoms_from_coord_file(xyz_path):
    """Infers the number of atoms from the first line of the coords file."""
    try:
        with open(xyz_path, "r") as f:
            first_line = f.readline()
            num_atoms = int(first_line.strip())
            logging.info(f"Inferred {num_atoms} atoms from first line of {xyz_path}.")
            return num_atoms
    except Exception as e:
        logging.warning(f"Could not infer atoms from {xyz_path}: {e}. Will try velocity file.")
        return None


# --- 3. Parallel Extraction Workers (Threaded) ---

def coords_extractor(args, want_frames, frame_map, chunk_size, stop_event, pbar, stats):
    """
    Worker thread to stream and extract coordinate frames.
    """
    logging.info("[CoordsWorker] Started.")
    try:
        for idx, (num_atoms, title, frame_data) in chunk_reader_coords(args.coords, chunk_size, stop_event):
            stats['coords_scanned'] += 1
            pbar.update(1)
            
            if idx in want_frames:
                try:
                    output_frame_num = frame_map[idx]['output_frame']
                    output_path = frame_map[idx]['xyz_path']
                    
                    with open(output_path, "w") as f:
                        f.write(f"{num_atoms}\n")
                        f.write(f"{title}\n") # Write original title line
                        for atom, x, y, z in frame_data:
                            f.write(f"{atom} {x:.10f} {y:.10f} {z:.10f}\n")
                    stats['coords_extracted'] += 1
                except KeyError:
                     # This can happen if frame_map isn't populated yet, but shouldn't
                     logging.warning(f"[CoordsWorker] Frame {idx} in want_frames but not in frame_map.")
                except Exception as e:
                    logging.error(f"[CoordsWorker] Failed to write {output_path}: {e}")
                    stats['errors'].append(f"Write failed: {output_path}")

            if idx > pbar.total:
                # Scanned past the expected last frame, update pbar total
                pbar.total = idx
                pbar.refresh()
            
            if stop_event.is_set():
                logging.info("[CoordsWorker] Stop event received.")
                break
                
    except Exception as e:
        logging.exception(f"[CoordsWorker] Fatal error: {e}")
        stats['errors'].append(f"CoordsWorker failed: {e}")
        stop_event.set() # Signal failure to other threads
    finally:
        pbar.close()
        logging.info(f"[CoordsWorker] Finished. Scanned: {stats['coords_scanned']}, Extracted: {stats['coords_extracted']}.")


def velocity_extractor(args, want_frames, frame_map, num_atoms, chunk_size, stop_event, pbar, stats):
    """
    Worker thread to stream and extract velocity frames.
    """
    if not args.velocities:
        logging.info("[VelWorker] No velocity file provided. Exiting.")
        pbar.close()
        return

    logging.info("[VelWorker] Started.")
    try:
        for idx, frame_data in chunk_reader_velocity(args.velocities, num_atoms, chunk_size, stop_event):
            stats['vel_scanned'] += 1
            pbar.update(1)
            
            if len(frame_data) != num_atoms:
                 logging.warning(f"[VelWorker] Frame {idx} has {len(frame_data)} atoms, expected {num_atoms}. Skipping.")
                 stats['mismatches'].append(f"Vel Frame {idx}: Atom count mismatch (found {len(frame_data)}, expected {num_atoms})")
                 continue
            
            if idx in want_frames:
                try:
                    output_frame_num = frame_map[idx]['output_frame']
                    output_path = frame_map[idx]['vel_path']

                    with open(output_path, "w") as f:
                        # Write in a simple .vel format (not XYZ)
                        f.write(f"{num_atoms}\n")
                        f.write(f"Source Frame {idx}\n")
                        for vx, vy, vz in frame_data:
                            f.write(f"{vx:.10f} {vy:.10f} {vz:.10f}\n")
                    stats['vel_extracted'] += 1
                except KeyError:
                     logging.warning(f"[VelWorker] Frame {idx} in want_frames but not in frame_map.")
                except Exception as e:
                    logging.error(f"[VelWorker] Failed to write {output_path}: {e}")
                    stats['errors'].append(f"Write failed: {output_path}")

            if idx > pbar.total:
                pbar.total = idx
                pbar.refresh()
                
            if stop_event.is_set():
                logging.info("[VelWorker] Stop event received.")
                break

    except Exception as e:
        logging.exception(f"[VelWorker] Fatal error: {e}")
        stats['errors'].append(f"VelWorker failed: {e}")
        stop_event.set() # Signal failure to other threads
    finally:
        pbar.close()
        logging.info(f"[VelWorker] Finished. Scanned: {stats['vel_scanned']}, Extracted: {stats['vel_extracted']}.")


# --- 4. Post-Processing Validation (Multiprocessing) ---

def get_specific_coord_frame(file_path, target_frame_index):
    """Reads a single, specific frame from the large coordinate file."""
    # Note: This is inefficient (streams from start) but is the only
    # robust way for non-indexed text files.
    chunk_size = 10 * 1024 * 1024 # 10MB chunks for validation streaming
    stop_event = Event() # Dummy event
    try:
        for idx, (num_atoms, title, frame_data) in chunk_reader_coords(file_path, chunk_size, stop_event):
            if idx == target_frame_index:
                return num_atoms, title, frame_data
        return None # Frame not found
    except Exception as e:
        logging.error(f"[Validation] Error reading coord frame {target_frame_index}: {e}")
        return None

def get_specific_velocity_frame(file_path, target_frame_index, num_atoms):
    """Reads a single, specific frame from the large velocity file."""
    chunk_size = 10 * 1024 * 1024 # 10MB chunks
    stop_event = Event() # Dummy event
    try:
        for idx, frame_data in chunk_reader_velocity(file_path, num_atoms, chunk_size, stop_event):
            if idx == target_frame_index:
                return frame_data
        return None # Frame not found
    except Exception as e:
        logging.error(f"[Validation] Error reading vel frame {target_frame_index}: {e}")
        return None

def read_extracted_xyz(file_path):
    """Reads a small, extracted .xyz file."""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            num_atoms = int(lines[0].strip())
            title = lines[1].strip() # .strip() is important
            data = []
            for line in lines[2:]:
                parts = line.split()
                if len(parts) == 4:
                    data.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
            return num_atoms, title, data
    except Exception as e:
        return f"Error: {e}"

def read_extracted_vel(file_path):
    """Reads a small, extracted .vel file."""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            num_atoms = int(lines[0].strip())
            data = []
            for line in lines[2:]:
                parts = line.split()
                if len(parts) == 3:
                    data.append(tuple(map(float, parts)))
            return num_atoms, data
    except Exception as e:
        return f"Error: {e}"

def _validate_frame_worker(args):
    """Internal worker for running validation on a single frame."""
    (
        output_frame_num,
        source_coord_idx,
        source_vel_idx,
        extracted_xyz_path,
        extracted_vel_path,
        source_coords_file,
        source_vel_file,
        num_atoms
    ) = args

    results = []

    # --- Validate Coords ---
    try:
        # 1. Read extracted file
        extracted_xyz = read_extracted_xyz(extracted_xyz_path)
        if isinstance(extracted_xyz, str):
             raise IOError(f"Failed to read extracted {extracted_xyz_path}: {extracted_xyz}")
        ext_n, ext_title, ext_data = extracted_xyz

        # 2. Read source file
        source_xyz = get_specific_coord_frame(source_coords_file, source_coord_idx)
        if source_xyz is None:
            raise IOError(f"Could not find source coord frame {source_coord_idx}")
        src_n, src_title, src_data = source_xyz
        
        # 3. Compare
        if ext_n != src_n:
            results.append(f"Frame {output_frame_num} [Coords]: Atom count mismatch (Ext: {ext_n}, Src: {src_n})")
        # We must .strip() the source title as it may have a newline
        if ext_title.strip() != src_title.strip():
             results.append(f"Frame {output_frame_num} [Coords]: Title line mismatch (Ext: '{ext_title}', Src: '{src_title}')")
        
        # Simple data comparison
        if len(ext_data) != len(src_data):
             results.append(f"Frame {output_frame_num} [Coords]: Data line count mismatch (Ext: {len(ext_data)}, Src: {len(src_data)})")
        elif str(ext_data) != str(src_data): # Quick and dirty comparison
             results.append(f"Frame {output_frame_num} [Coords]: Data mismatch (content differs)")
        
    except Exception as e:
        results.append(f"Frame {output_frame_num} [Coords] VALIDATION FAILED: {e}")

    # --- Validate Velocities ---
    if source_vel_file:
        try:
            # 1. Read extracted file
            extracted_vel = read_extracted_vel(extracted_vel_path)
            if isinstance(extracted_vel, str):
                raise IOError(f"Failed to read extracted {extracted_vel_path}: {extracted_vel}")
            ext_n_vel, ext_data_vel = extracted_vel

            # 2. Read source file
            source_vel_data = get_specific_velocity_frame(source_vel_file, source_vel_idx, num_atoms)
            if source_vel_data is None:
                raise IOError(f"Could not find source vel frame {source_vel_idx}")
            src_n_vel = len(source_vel_data)

            # 3. Compare
            if ext_n_vel != src_n_vel:
                results.append(f"Frame {output_frame_num} [Vels]: Atom count mismatch (Ext: {ext_n_vel}, Src: {src_n_vel})")
            elif str(ext_data_vel) != str(source_vel_data):
                results.append(f"Frame {output_frame_num} [Vels]: Data mismatch (content differs)")

        except Exception as e:
            results.append(f"Frame {output_frame_num} [Vels] VALIDATION FAILED: {e}")
    
    if not results:
        return (output_frame_num, "OK")
    else:
        return (output_frame_num, "MISMATCH", results)


def run_validation(args, frame_map, num_atoms, num_workers, report_file):
    """
    Manages the post-processing validation.
    """
    logging.info("--- Starting Post-Processing Validation ---")
    
    frames_to_validate = []
    all_output_frames = sorted(frame_map.keys(), key=lambda x: frame_map[x]['output_frame'])
    
    if args.full_validation:
        frames_to_validate = all_output_frames
        logging.info(f"Running FULL validation on {len(frames_to_validate)} extracted frames.")
    else:
        sample_size = min(len(all_output_frames), VALIDATION_SAMPLE_SIZE)
        frames_to_validate = random.sample(all_output_frames, sample_size)
        logging.info(f"Running RANDOM validation on {sample_size} extracted frames.")

    if not frames_to_validate:
        logging.warning("No frames were extracted. Skipping validation.")
        return 0, 0

    # Prepare arguments for pool.map
    validation_tasks = []
    for source_idx in frames_to_validate:
        frame_info = frame_map[source_idx]
        validation_tasks.append((
            frame_info['output_frame'],
            frame_info['source_coord_idx'],
            frame_info['source_vel_idx'],
            frame_info['xyz_path'],
            frame_info['vel_path'],
            args.coords,
            args.velocities,
            num_atoms
        ))

    mismatches = 0
    validated_ok = 0
    
    with open(report_file, 'w') as f:
        f.write("--- Validation Report ---\n")
        f.write(f"Timestamp: {datetime.now().isoformat()}\n")
        f.write(f"Mode: {'Full' if args.full_validation else f'Random Sample ({len(frames_to_validate)})'}\n\n")
        
        logging.info(f"Starting validation with {num_workers} workers...")
        with tqdm(total=len(validation_tasks), desc="Validating Frames", unit="frame") as pbar:
            with multiprocessing.Pool(processes=num_workers) as pool:
                for result in pool.imap_unordered(_validate_frame_worker, validation_tasks):
                    pbar.update(1)
                    if result[1] == "OK":
                        validated_ok += 1
                    else:
                        mismatches += 1
                        f.write(f"--- MISMATCH: Frame {result[0]} (Source: {frame_map[frames_to_validate[pbar.n-1]]['source_coord_idx']}) ---\n")
                        for line in result[2]:
                            f.write(f"  - {line}\n")
                        f.write("\n")

    logging.info(f"Validation complete. OK: {validated_ok}, Mismatches: {mismatches}.")
    logging.info(f"Validation report saved to {report_file}")
    
    return validated_ok, mismatches


# --- 5. Resource Management & Setup ---

def get_system_resources(args_workers):
    """Gets available system resources and determines worker counts."""
    cpu_cores = multiprocessing.cpu_count()
    available_mem = psutil.virtual_memory().available
    available_mem_gb = available_mem / (1024**3)
    
    # Determine workers for counting/validation
    # Leave one core for system/main thread
    validation_workers = max(1, cpu_cores - 1)
    if args_workers:
        validation_workers = min(args_workers, cpu_cores)
        
    logging.info(f"System Resources: {cpu_cores} CPU Cores, {available_mem_gb:.2f} GB RAM available.")
    logging.info(f"Will use {validation_workers} workers for counting and validation.")
    
    if available_mem_gb < 2.0:
         logging.warning(f"Low available memory ({available_mem_gb:.2f} GB). Processing large files may fail.")
         
    return validation_workers, available_mem

def setup_logging(output_dir, log_level_str):
    """Configures logging to file and console."""
    
    log_level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
    }
    log_level = log_level_map.get(log_level_str.upper(), logging.INFO)

    # Ensure output dir exists
    os.makedirs(output_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(output_dir, f"processing_log_{timestamp}.log")
    
    # Get root logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG) # Capture all levels
    
    # Remove any existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Console handler (respects user's level)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_format = logging.Formatter("[%(levelname)s] %(message)s")
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)

    # File handler (always DEBUG)
    file_handler = logging.FileHandler(log_file, 'w')
    file_handler.setLevel(logging.DEBUG)
    file_format = logging.Formatter("%(asctime)s [%(levelname)s] (%(threadName)s) %(message)s")
    file_handler.setFormatter(file_format)
    logger.addHandler(file_handler)
    
    logging.info(f"Logging configured. Console level: {log_level_str.upper()}, Log file: {log_file}")
    return timestamp
    

# --- 6. Main Execution ---

def main():
    parser = argparse.ArgumentParser(
        description="Production-Ready MD Trajectory Processor",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # --- Input/Output Args ---
    parser.add_argument("-xyz", "--coords", type=str, required=True, help="Path to the coords.xyz file.")
    parser.add_argument("-vel", "--velocities", type=str, help="Path to the velocity.out file (optional).")
    parser.add_argument("-o", "--output", type=str, default="extracted_frames", help="Base output directory (default: extracted_frames).")

    # --- Frame Selection Args ---
    parser.add_argument("--start", type=int, default=1, help="Starting frame index (1-based). Default: 1.")
    parser.add_argument("--stop", type=int, help="Stopping frame index (inclusive). Default: all frames.")
    parser.add_argument("--step", type=int, default=1, help="Step size for frame extraction (default: 1).")

    # --- Performance & Validation Args ---
    parser.add_argument("--workers", type=int, help="Max CPU cores for counting/validation (default: all-1).")
    parser.add_argument("--chunk-size-mb", type=int, default=DEFAULT_CHUNK_MB, help=f"Chunk size in MB for reading (default: {DEFAULT_CHUNK_MB}).")
    parser.add_argument("--full-validation", action="store_true", help=f"Validate ALL extracted frames, not just a sample of {VALIDATION_SAMPLE_SIZE}.")
    parser.add_argument("--skip-validation", action="store_true", help="Skip the entire post-processing validation step.")
    
    # --- Utility Args ---
    parser.add_argument("--dry-run", action="store_true", help="Count frames, plan, and report without extracting files.")
    parser.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Set console logging level (default: INFO).")

    args = parser.parse_args()

    # --- 1. Setup Logging and Environment ---
    start_time = time.time()
    timestamp = setup_logging(args.output, args.log_level)
    print_banner()
    
    logging.info(f"Script started at {datetime.now().isoformat()}")
    logging.info(f"Run arguments: {vars(args)}")

    # Shared stop event for graceful shutdown
    stop_event = Event()
    
    # Statistics dictionary
    stats = {
        'coords_scanned': 0,
        'coords_extracted': 0,
        'vel_scanned': 0,
        'vel_extracted': 0,
        'errors': [],
        'mismatches': [],
        'validation_ok': 0,
        'validation_mismatch': 0,
    }
    
    # Summary report file
    summary_file = os.path.join(args.output, f"summary_report_{timestamp}.txt")
    validation_report_file = os.path.join(args.output, f"validation_report_{timestamp}.txt")
    
    # Main try/finally block to ensure summary is always written
    try:
        # --- 2. Resource Assessment ---
        num_workers, available_mem = get_system_resources(args.workers)
        
        # --- 3. Parameter Validation and Planning ---
        if not os.path.isfile(args.coords):
            raise FileNotFoundError(f"Coords file not found: {args.coords}")
        if args.velocities and not os.path.isfile(args.velocities):
             raise FileNotFoundError(f"Velocity file not found: {args.velocities}")

        if args.start < 1:
            raise ValueError("--start must be >= 1 (1-based indexing).")
        if args.step < 1:
            raise ValueError("--step must be >= 1.")

        # --- 4. Count Total Frames ---
        # This is the primary requirement for robust frame detection
        logging.info(f"Counting total frames in {args.coords} using pattern '{FRAME_MARKER}'...")
        n_frames = count_pattern_in_file(args.coords, FRAME_MARKER, num_workers)
        
        if n_frames == 0:
            # Fallback: Try to infer from coord file structure (slower)
            logging.warning(f"No occurrences of '{FRAME_MARKER}' found.")
            logging.warning("Attempting to count frames by parsing... This may be slow and/or inaccurate.")
            n_frames_parsed = 0
            # Use a dummy pbar
            dummy_pbar = tqdm(total=0, desc="Parsing Frames", unit="frame", leave=False)
            for _ in chunk_reader_coords(args.coords, args.chunk_size_mb * 1024 * 1024, stop_event):
                n_frames_parsed += 1
                dummy_pbar.total = n_frames_parsed
                dummy_pbar.refresh()
            dummy_pbar.close()
            
            if n_frames_parsed > 0:
                logging.warning(f"Found {n_frames_parsed} frames by parsing. Using this as total.")
                n_frames = n_frames_parsed
            else:
                 raise ValueError("Could not find any frames in the coordinate file.")
        
        # --- 5. Finalize Frame Logic and Resource Scaling ---
        start_frame = args.start
        step = args.step
        stop_frame = args.stop if args.stop is not None else n_frames

        if args.stop is not None and args.stop > n_frames:
            logging.warning(f"--stop ({args.stop}) is greater than detected frames ({n_frames}). Capping to {n_frames}.")
            stop_frame = n_frames
            
        if stop_frame < start_frame:
            raise ValueError(f"--stop ({stop_frame}) must be >= --start ({start_frame}).")
            
        # Build the set of frames we want to extract
        want_frames_list = list(range(start_frame, stop_frame + 1, step))
        want_frames_set = set(want_frames_list)
        total_to_extract = len(want_frames_list)
        
        logging.info(f"Total frames in file: {n_frames}")
        logging.info(f"Extraction plan: Start={start_frame}, Stop={stop_frame}, Step={step}")
        logging.info(f"Total frames to extract: {total_to_extract}")

        # Dynamic Resource Management
        chunk_size_bytes = args.chunk_size_mb * 1024 * 1024
        if n_frames >= LARGE_FILE_FRAMES:
            # Increase chunk size for very large files if we have memory
            if available_mem > 4 * 1024**3: # > 4GB RAM
                new_chunk_mb = min(512, args.chunk_size_mb * 4) # Cap at 512MB
                chunk_size_bytes = new_chunk_mb * 1024 * 1024
                logging.info(f"Large trajectory detected. Increasing chunk size to {new_chunk_mb} MB.")
            else:
                 logging.warning("Large trajectory detected but available memory is low. Sticking with {args.chunk_size_mb} MB chunk size.")

        # --- 6. Generate Frame Map and Output Dirs ---
        logging.info("Generating frame map and creating directories...")
        frame_map_file = os.path.join(args.output, "coords_frame_indices_map.csv")
        frame_map = {} # dict for fast lookup: {source_idx: info_dict}
        
        with open(frame_map_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["OutputFrame", "SourceCoordsIndex", "SourceVelocityIndex"])
            
            pbar_map = tqdm(total=total_to_extract, desc="Preparing Output", unit="frame", leave=False)
            for i, source_frame_index in enumerate(want_frames_list):
                output_frame_num = i + 1

                # --- New Directory Structure ---
                # Subdirectory name: frame_1, frame_2, etc.
                subdir_name = f"frame_{output_frame_num}"
                output_subdir = os.path.join(args.output, subdir_name)
                
                # Create the subdirectory if not in dry run
                if not args.dry_run:
                    os.makedirs(output_subdir, exist_ok=True)
                
                # File name: frame_1.xyz, frame_2.xyz, etc.
                output_file_name_base = f"frame_{output_frame_num}"
                output_xyz = os.path.join(output_subdir, f"{output_file_name_base}.xyz")
                output_vel = os.path.join(output_subdir, f"{output_file_name_base}.vel")
                
                # Write to CSV (use the base name for the map)
                writer.writerow([output_file_name_base, source_frame_index, source_frame_index])
                
                # Add to lookup map
                frame_map[source_frame_index] = {
                    'output_frame': output_frame_num,
                    'source_coord_idx': source_frame_index,
                    'source_vel_idx': source_frame_index, # Assume 1:1 mapping
                    'xyz_path': output_xyz,
                    'vel_path': output_vel,
                }
                pbar_map.update(1)
            pbar_map.close()
        
        logging.info(f"Frame map saved to {frame_map_file}")
        
        # --- 7. Dry Run ---
        if args.dry_run:
            logging.info("--- DRY RUN ENABLED ---")
            logging.info("Processing plan complete. No files will be extracted.")
            logging.info(f"Output directory: {args.output}")
            logging.info(f"Frames to extract: {total_to_extract}")
            logging.info(f"Log file: {os.path.join(args.output, f'processing_log_{timestamp}.log')}")
            logging.info(f"Frame map: {frame_map_file}")
            logging.info("Dry run complete. Exiting.")
            # Set to 0 so summary report is accurate
            stats['coords_scanned'] = n_frames 
            stats['vel_scanned'] = n_frames if args.velocities else 0
            return # Exit before extraction

        # --- 8. Infer Atoms per Frame ---
        # This is critical for the velocity parser
        num_atoms = infer_atoms_from_coord_file(args.coords)
        if num_atoms is None and args.velocities:
            num_atoms = infer_atoms_from_velocity_file(args.velocities)
            
        if num_atoms is None:
            # This is only a problem if velocities are requested
            if args.velocities:
                 raise ValueError("Could not determine number of atoms per frame from coords or velocity file. Cannot process velocities.")
            else:
                logging.warning("Could not determine number of atoms. Proceeding with coords-only extraction.")
                # The coord parser will determine atoms dynamically
        
        # --- 9. Start Parallel Extraction ---
        logging.info("--- Starting Parallel Extraction ---")
        
        # TQDM progress bars
        # Note: We set total to n_frames, as the scanners
        # will read the *entire* file.
        coords_pbar = tqdm(total=n_frames, desc="Coords Processed", unit="frame", position=0)
        vel_pbar = tqdm(total=n_frames if args.velocities else 0, desc="Vels Processed  ", unit="frame", position=1)
        
        # Start worker threads
        coords_thread = Thread(
            target=coords_extractor,
            name="CoordsWorker",
            args=(args, want_frames_set, frame_map, chunk_size_bytes, stop_event, coords_pbar, stats),
            daemon=True
        )
        
        vel_thread = Thread(
            target=velocity_extractor,
            name="VelWorker",
            args=(args, want_frames_set, frame_map, num_atoms, chunk_size_bytes, stop_event, vel_pbar, stats),
            daemon=True
        )
        
        coords_thread.start()
        vel_thread.start()
        
        # Wait for threads to finish, handling KeyboardInterrupt
        while coords_thread.is_alive() or vel_thread.is_alive():
            time.sleep(0.5)

        logging.info("Extraction threads have finished.")

        # --- 10. Final Scan Check ---
        if args.velocities:
            if stats['coords_scanned'] != stats['vel_scanned']:
                msg = (f"Frame scan count mismatch! "
                       f"Coords: {stats['coords_scanned']}, Vels: {stats['vel_scanned']}")
                logging.warning(msg)
                stats['mismatches'].append(msg)
            
        if stats['coords_extracted'] != total_to_extract:
             msg = (f"Coordinate extraction mismatch! "
                    f"Expected: {total_to_extract}, Extracted: {stats['coords_extracted']}")
             logging.warning(msg)
             stats['mismatches'].append(msg)

        if args.velocities and stats['vel_extracted'] != total_to_extract:
             msg = (f"Velocity extraction mismatch! "
                    f"Expected: {total_to_extract}, Extracted: {stats['vel_extracted']}")
             logging.warning(msg)
             stats['mismatches'].append(msg)
             
        # --- 11. Run Validation ---
        if args.skip_validation:
            logging.info("Skipping validation as requested.")
        elif num_atoms is None:
            logging.error("Skipping validation: could not determine num_atoms.")
            stats['errors'].append("Validation skipped (unknown num_atoms)")
        elif stats['coords_extracted'] == 0:
            logging.warning("Skipping validation: no coordinate frames were extracted.")
        else:
            ok, mismatch = run_validation(
                args, frame_map, num_atoms, num_workers, validation_report_file
            )
            stats['validation_ok'] = ok
            stats['validation_mismatch'] = mismatch
            if mismatch > 0:
                 stats['errors'].append(f"{mismatch} validation mismatches found.")

    except (KeyboardInterrupt, SystemExit):
        logging.warning("\n--- INTERRUPT RECEIVED ---")
        logging.warning("Signaling threads to stop. Writing partial summary...")
        stop_event.set() # Tell threads to stop
        stats['errors'].append("User Interrupted")
    except FileNotFoundError as e:
        logging.error(f"File Not Found: {e}")
        stats['errors'].append(str(e))
    except ValueError as e:
        logging.error(f"Value Error: {e}")
        stats['errors'].append(str(e))
    except Exception as e:
        logging.exception(f"An unexpected error occurred: {e}")
        stop_event.set()
        stats['errors'].append(f"Fatal: {e}")
    finally:
        # --- 12. Write Summary Report ---
        elapsed = time.time() - start_time
        logging.info("--- Processing Summary ---")
        logging.info(f"Total processing time: {elapsed:.2f} seconds")
        
        proc = psutil.Process(os.getpid())
        peak_mem = proc.memory_info().rss / (1024**2) # in MB
        
        summary_lines = [
            "========================================",
            "        Processing Summary Report",
            "========================================",
            f"Timestamp: {datetime.now().isoformat()}",
            f"Total Time: {elapsed:.2f} seconds",
            f"Peak Memory Usage: {peak_mem:.2f} MB",
            "",
            "--- Configuration ---",
            f"Coords File: {args.coords}",
            f"Velocity File: {args.velocities}",
            f"Output Directory: {args.output}",
            f"Plan: Start={start_frame}, Stop={stop_frame}, Step={step}",
            f"Dry Run: {args.dry_run}",
            "",
            "--- Extraction Statistics ---",
            f"Total Frames in Coords File: {stats.get('coords_scanned', 'N/A')}",
            f"Total Frames in Velocity File: {stats.get('vel_scanned', 'N/A')}",
            f"Frames Requested to Extract: {total_to_extract if not args.dry_run else 'N/A (Dry Run)'}",
            f"Frames Successfully Extracted (Coords): {stats['coords_extracted']}",
            f"Frames Successfully Extracted (Vels): {stats['vel_extracted']}",
            "",
            "--- Validation Statistics ---",
            f"Validation Skipped: {args.skip_validation}",
            f"Validation Mode: {'Full' if args.full_validation else f'Random Sample (Max {VALIDATION_SAMPLE_SIZE})'}",
            f"Frames Validated OK: {stats['validation_ok']}",
            f"Frames with Mismatches: {stats['validation_mismatch']}",
            f"Validation Report: {validation_report_file if not args.skip_validation else 'N/A'}",
            "",
            "--- Errors and Mismatches ---",
        ]
        
        if not stats['errors'] and not stats['mismatches']:
            summary_lines.append("No errors or mismatches detected.")
            logging.info("Processing complete. No errors.")
        else:
            logging.error(f"Processing finished with {len(stats['errors'])} errors and {len(stats['mismatches'])} mismatches.")
            if stats['mismatches']:
                summary_lines.append("Mismatch warnings (see log for details):")
                for msg in stats['mismatches']:
                    summary_lines.append(f"  - {msg}")
            if stats['errors']:
                summary_lines.append("Errors (see log for details):")
                for err in stats['errors']:
                    summary_lines.append(f"  - {err}")

        summary_lines.append("========================================")
        
        try:
            with open(summary_file, 'w') as f:
                f.write("\n".join(summary_lines))
            logging.info(f"Summary report saved to {summary_file}")
        except Exception as e:
            logging.error(f"Failed to write summary report: {e}")
            
        print(f"\nProcessing complete. See log and summary files in {args.output}")

if __name__ == "__main__":
    # Set multiprocessing start method for compatibility (esp. on macOS/Windows)
    try:
        multiprocessing.set_start_method("spawn")
    except RuntimeError:
        pass # Already set
        
    main()


