import os
import sys
import argparse
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
import logging
import time
from tqdm import tqdm
from threading import Thread


# Banner for the script
def print_banner():
    logging.info("=========================")
    logging.info("Trajectory Processor Script")
    logging.info("=========================")


# Required Functions
def process_velocity_chunk(chunk, leftover="", atoms_per_frame=0):
    results = []
    lines = leftover + chunk  # Include leftover from the previous chunk
    lines = lines.splitlines()
    frame = []
    in_velocity_block = False

    for line in lines:
        # Identify the start and end of a velocity block
        if "$VELOC" in line:
            in_velocity_block = True
            logging.debug("process_velocity_chunk: enter $VELOC block")
            continue
        if "$ENDVELOC" in line:
            in_velocity_block = False
            logging.debug(
                f"process_velocity_chunk: exit $VELOC block; current frame length={len(frame)}"
            )
            # flush frame if it matches expected count or if atoms_per_frame not specified
            if atoms_per_frame == 0 or len(frame) == atoms_per_frame:
                if frame:
                    results.append(frame)
                    frame = []
            continue

        # Skip lines outside velocity blocks
        if not in_velocity_block:
            continue

        # Parse velocity values
        parts = line.split()
        if len(parts) == 3:
            try:
                frame.append(tuple(map(float, parts)))
            except ValueError:
                logging.debug(f"Skipping malformed velocity line: {line}")

        # If the current frame is complete (and atoms_per_frame specified), append it
        if atoms_per_frame and len(frame) == atoms_per_frame:
            logging.debug(
                f"process_velocity_chunk: complete frame of length {len(frame)}"
            )
            results.append(frame)
            frame = []

    # Return leftover lines that could not form a complete frame (store numeric text)
    leftover = "\n".join([f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}" for v in frame])
    return results, leftover


def chunk_reader_velocity(
    file_path, chunk_size=1024 * 1024 * 50, atoms_per_frame=0
):  # 50MB per chunk
    with open(file_path, "r") as file:
        leftover = ""
        while True:
            chunk = file.read(chunk_size)
            if not chunk:
                logging.debug(
                    f"chunk_reader_velocity: EOF reached, leftover length={len(leftover)}"
                )
                # no more data; attempt to flush leftover (process any remaining complete frames)
                if leftover:
                    results, leftover = process_velocity_chunk(
                        "", leftover, atoms_per_frame
                    )
                    if results:
                        logging.debug(
                            f"chunk_reader_velocity: yielding {len(results)} frames from leftover"
                        )
                        yield results, leftover
                break
            # debug preview of chunk to ensure markers are present
            logging.debug(f"chunk_reader_velocity: chunk preview: {repr(chunk[:200])}")
            results, leftover = process_velocity_chunk(chunk, leftover, atoms_per_frame)
            logging.debug(
                f"chunk_reader_velocity: read chunk size={len(chunk)}, yielded {len(results)} frames, leftover length={len(leftover)}"
            )
            yield results, leftover


def process_velocity_file_in_parallel(file_path, atoms_per_frame, num_processes=8):
    # For very large velocity files it's safer to stream; this helper retains previous behavior
    pool = multiprocessing.Pool(num_processes)
    chunks_with_leftovers = list(
        chunk_reader_velocity(file_path, atoms_per_frame=atoms_per_frame)
    )
    all_results = []
    for results, _ in chunks_with_leftovers:
        all_results.extend(results)
    pool.close()
    pool.join()
    return all_results


def extract_frames_as_velocities(
    velocity_data, start=0, stop=None, step=1, output_dir="extracted_frames"
):
    """
    Extracts specific frames of velocities and saves them as individual .vel files.

    Parameters:
        velocity_data (list): List of frames where each frame is a list of velocity tuples.
        start (int): Starting frame index (1-based).
        stop (int): Stopping frame index (inclusive).
        step (int): Step size for extraction (default is 1).
        output_dir (str): Directory to save the extracted velocities (default is 'extracted_frames').
    """
    # Set default stop value if not provided
    if stop is None:
        stop = len(velocity_data)

    # Validate inputs
    if start < 1 or stop > len(velocity_data) or step <= 0:
        raise ValueError("Invalid start, stop, or step values for frame extraction.")

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Extract velocity frames based on user-specified range
    frame_indices = range(start, stop + 1, step)
    for idx, frame_num in enumerate(frame_indices):
        os.makedirs(f"{output_dir}/frame_{idx + 1}", exist_ok=True)
        frame = velocity_data[frame_num - 1]
        file_path = os.path.join(
            f"{output_dir}/frame_{idx + 1}", f"frame_{idx + 1}.vel"
        )
        with open(file_path, "w") as f:
            f.write(f"{len(frame)}\n")
            f.write(f"Frame {frame_num}\n")
            for v in frame:
                f.write(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n")

    return output_dir


def extract_frames_as_xyz(
    trajectory_data, start=0, stop=None, step=1, output_dir="extracted_frames"
):
    """
    Extracts specific frames from the trajectory and saves them as individual .xyz files.

    Parameters:
        trajectory_data (list): List of frames where each frame is a list of atom data.
        start (int): Starting frame index (1-based).
        stop (int): Stopping frame index (inclusive).
        step (int): Step size for extraction (default is 1).
        output_dir (str): Directory to save the extracted frames (default is 'extracted_frames').
    """
    # Set default stop value if not provided
    if stop is None:
        stop = len(trajectory_data)

    # Validate inputs
    if start < 1 or stop > len(trajectory_data) or step <= 0:
        raise ValueError("Invalid start, stop, or step values for frame extraction.")

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Extract frames based on user-specified range
    frame_indices = range(start, stop + 1, step)
    for idx, frame_num in enumerate(frame_indices):
        os.makedirs(f"{output_dir}/frame_{idx + 1}", exist_ok=True)
        frame = trajectory_data[frame_num - 1]
        file_path = os.path.join(
            f"{output_dir}/frame_{idx + 1}", f"frame_{idx + 1}.xyz"
        )
        logging.debug(f"Writing {file_path}")
        with open(file_path, "w") as f:
            f.write(f"{len(frame)}\n")
            f.write(f"Frame {frame_num}\n")
            for atom in frame:
                f.write(f"{atom[0]} {atom[1]:.10f} {atom[2]:.10f} {atom[3]:.10f}\n")

    return output_dir


# Function to process a chunk of the trajectory file
def process_chunk(chunk, leftover=""):
    results = []
    lines = leftover + chunk  # Include leftover from the previous chunk
    lines = lines.splitlines()
    frame_size = 0
    i = 0

    while i < len(lines):
        try:
            # Dynamically determine the number of atoms (frame size)
            num_atoms = int(lines[i].strip())
            frame_size = num_atoms + 2  # Include the header and metadata
        except ValueError:
            i += 1
            continue

        if i + frame_size > len(lines):  # Incomplete frame
            leftover = "\n".join(lines[i:])
            break

        frame_data = lines[i : i + frame_size]
        frame_results = []
        for atom_line in frame_data[2:]:  # Skip the first two lines (header + metadata)
            parts = atom_line.split()
            if len(parts) == 4:  # Ensure line has atomic symbol and 3 coordinates
                try:
                    atom_symbol, x, y, z = parts
                    frame_results.append((atom_symbol, float(x), float(y), float(z)))
                except ValueError:
                    logging.debug(f"Skipping malformed atom line: {atom_line}")

        results.append(frame_results)
        i += frame_size

    return results, leftover


# Function to split the file into manageable chunks for multiprocessing
def chunk_reader(file_path, chunk_size=1024 * 1024 * 50):  # 50MB per chunk
    with open(file_path, "r") as file:
        leftover = ""
        while True:
            chunk = file.read(chunk_size)
            if not chunk:
                # no more data; attempt to flush leftover (process any remaining complete frames)
                if leftover:
                    results, leftover = process_chunk("", leftover)
                    if results:
                        yield results, leftover
                break
            results, leftover = process_chunk(chunk, leftover)
            yield results, leftover


# Multiprocessing function
def process_file_in_parallel(file_path, num_processes=8):
    pool = multiprocessing.Pool(num_processes)
    chunks_with_leftovers = list(chunk_reader(file_path))
    all_results = []
    for results, _ in chunks_with_leftovers:
        all_results.extend(results)
    pool.close()
    pool.join()
    return all_results


# --- New utilities for streaming extraction and counting ---


def _count_worker(args):
    """Top-level worker for counting pattern occurrences in a byte range.

    args: tuple(file_path, pattern, start, end)
    """
    file_path, pattern, start, end = args
    overlap = len(pattern) - 1
    read_start = max(0, start - overlap)
    read_end = end + overlap
    read_len = read_end - read_start
    count = 0
    with open(file_path, "r", errors="ignore") as f:
        f.seek(read_start)
        data = f.read(read_len)
        idx = data.find(pattern)
        while idx != -1:
            abs_pos = read_start + idx
            if abs_pos >= start and abs_pos < end:
                count += 1
            idx = data.find(pattern, idx + 1)
    return count


def count_pattern_in_file(
    file_path, pattern="FINAL HEAT OF FORMATION =", num_workers=None
):
    """
    Count occurrences of a pattern in a large text file using multiple processes.
    Returns integer count.
    """
    if num_workers is None:
        num_workers = max(1, multiprocessing.cpu_count() - 1)

    file_size = os.path.getsize(file_path)
    if file_size == 0:
        return 0

    pattern_len = len(pattern)
    part_size = max(
        10 * 1024 * 1024, file_size // num_workers + 1
    )  # at least 10MB per part
    ranges = []
    for i in range(0, file_size, part_size):
        start = i
        end = min(file_size, i + part_size)
        ranges.append((start, end))

    # Prepare args for top-level worker to avoid pickling local functions
    worker_args = [(file_path, pattern, start, end) for (start, end) in ranges]
    with multiprocessing.Pool(min(num_workers, len(ranges))) as pool:
        counts = pool.map(_count_worker, worker_args)

    total = sum(counts)
    return total


def stream_coords_frames(coords_path, chunk_size=1024 * 1024 * 50):
    """Generator yielding (frame_index, frame_data) for coordinates file. frame_index starts at 1."""
    frame_idx = 0
    for results, _ in chunk_reader(coords_path, chunk_size=chunk_size):
        for frame in results:
            frame_idx += 1
            yield frame_idx, frame


def stream_velocity_frames(vel_path, atoms_per_frame, chunk_size=1024 * 1024 * 50):
    """Generator yielding (frame_index, frame_data) for velocities file. frame_index starts at 1."""
    frame_idx = 0
    for results, _ in chunk_reader_velocity(
        vel_path, chunk_size=chunk_size, atoms_per_frame=atoms_per_frame
    ):
        for frame in results:
            frame_idx += 1
            yield frame_idx, frame


def infer_atoms_from_velocity_file(vel_path, max_lines=100000):
    """Attempt to infer the number of atoms per velocity frame by scanning for first $VELOC block."""
    in_block = False
    count = 0
    try:
        with open(vel_path, "r") as f:
            for i, line in enumerate(f):
                if i > max_lines:
                    break
                if "$VELOC" in line:
                    in_block = True
                    count = 0
                    continue
                if "$ENDVELOC" in line and in_block:
                    return count
                if in_block:
                    parts = line.split()
                    if len(parts) == 3:
                        count += 1
    except Exception:
        return None
    return None


# Main script
if __name__ == "__main__":
    # Basic console logger until output dir is known. Set to DEBUG to capture parser debug info.
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s [%(levelname)s] %(message)s"
    )

    print_banner()

    parser = argparse.ArgumentParser(description="Process molecular trajectory files.")
    parser.add_argument(
        "-xyz", "--coords", type=str, required=True, help="Path to the coords.xyz file."
    )
    parser.add_argument(
        "-vel", "--velocities", type=str, help="Path to the velocity.out file."
    )
    parser.add_argument(
        "--start",
        type=int,
        default=1,
        help="Starting frame index (1-based). Default: 1.",
    )
    parser.add_argument(
        "--stop",
        type=int,
        help="Stopping frame index (inclusive). Default: computed from coords file.",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=1,
        help="Step size for frame extraction (default: 1).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="extracted_frames",
        help="Base output directory (default: extracted_frames).",
    )
    parser.add_argument(
        "--count-workers",
        type=int,
        help="Number of workers to use when counting frames.",
    )
    parser.add_argument(
        "--chunk-size-mb",
        type=int,
        default=50,
        help="Chunk size in MB for reading (default:50).",
    )

    args = parser.parse_args()

    start_time = time.time()
    summary = {
        "coords_scanned": 0,
        "coords_extracted": 0,
        "vel_scanned": 0,
        "vel_extracted": 0,
        "errors": [],
    }

    try:
        if not os.path.isfile(args.coords):
            raise FileNotFoundError(f"Coords file not found: {args.coords}")

        # Read first line to get number of atoms per frame (if standard formatting)
        try:
            with open(args.coords, "r") as cf:
                first_line = cf.readline()
                num_atoms_per_frame = int(first_line.strip())
        except Exception as e:
            logging.warning(f"Could not read atoms-per-frame from first line: {e}")
            num_atoms_per_frame = None

        # Determine number of frames (stop) if not provided
        needed_to_count = args.stop is None
        if needed_to_count:
            logging.info(
                "Counting frames in coords file by searching for 'FINAL HEAT OF FORMATION =' occurrences..."
            )
            n_frames = count_pattern_in_file(
                args.coords,
                pattern="FINAL HEAT OF FORMATION =",
                num_workers=args.count_workers,
            )
            logging.info(f"Found {n_frames} frames (pattern occurrences).")
            if n_frames == 0:
                # Fallback: try to estimate from streaming parse
                logging.info(
                    "No pattern occurrences found; attempting to estimate frames by streaming parsing (may be slow)..."
                )
                n_frames = 0
                for _ in stream_coords_frames(
                    args.coords, chunk_size=args.chunk_size_mb * 1024 * 1024
                ):
                    n_frames += 1
                logging.info(f"Estimated {n_frames} frames from parsing.")
        else:
            n_frames = None

        # Resolve start/stop/step
        start = args.start if args.start is not None else 1
        step = args.step if args.step is not None else 1
        stop = args.stop if args.stop is not None else n_frames

        if stop is None:
            raise ValueError(
                "Could not determine stop frame. Provide --stop or ensure the coords file contains the pattern."
            )

        if start < 1:
            raise ValueError("--start must be >= 1 (1-based indexing).")
        if stop < start:
            raise ValueError("--stop must be >= --start.")

        # Cap stop to n_frames if known
        if n_frames is not None and stop > n_frames:
            logging.warning(
                f"--stop ({stop}) is greater than detected frames ({n_frames}). Capping to {n_frames}."
            )
            stop = n_frames

        total_frames = n_frames if n_frames is not None else stop
        logging.info(
            f"Extraction plan: start={start}, stop={stop}, step={step}, total_frames={total_frames}"
        )

        # Prepare output dir and logging to file
        os.makedirs(args.output, exist_ok=True)
        fh = logging.FileHandler(os.path.join(args.output, "processing.log"))
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logging.getLogger().addHandler(fh)

        # Build set for quick membership check
        want_frames = set(range(start, stop + 1, step))
        # Precompute mapping from frame number -> extraction index to avoid repeated sorted/index calls
        sorted_want = sorted(want_frames)
        extraction_map = {frame: i + 1 for i, frame in enumerate(sorted_want)}

        # Threads to stream and write frames
        coords_pbar = tqdm(
            total=total_frames, desc="Coords processed", unit="frame", position=0
        )
        vel_pbar = tqdm(
            total=total_frames if args.velocities else 0,
            desc="Velocities processed",
            unit="frame",
            position=1,
        )

        def coords_worker():
            try:
                for idx, frame in stream_coords_frames(
                    args.coords, chunk_size=args.chunk_size_mb * 1024 * 1024
                ):
                    summary["coords_scanned"] += 1
                    coords_pbar.update(1)
                    if idx in want_frames:
                        extraction_index = extraction_map[idx]
                        out_sub = os.path.join(args.output, f"frame_{extraction_index}")
                        os.makedirs(out_sub, exist_ok=True)
                        file_path = os.path.join(
                            out_sub, f"frame_{extraction_index}.xyz"
                        )
                        with open(file_path, "w") as f:
                            f.write(f"{len(frame)}\n")
                            f.write(f"Frame {idx}\n")
                            for atom in frame:
                                f.write(
                                    f"{atom[0]} {atom[1]:.10f} {atom[2]:.10f} {atom[3]:.10f}\n"
                                )
                        summary["coords_extracted"] += 1
            except Exception as e:
                logging.exception(f"Error while streaming coords: {e}")
                summary["errors"].append(str(e))
            finally:
                coords_pbar.close()

        def vel_worker():
            if not args.velocities:
                vel_pbar.close()
                return
            if num_atoms_per_frame is None:
                # try to infer a reasonable default
                logging.warning(
                    "Number of atoms per frame unknown; attempting to read first velocity block to infer atoms per frame."
                )
                inferred = infer_atoms_from_velocity_file(args.velocities)
                if inferred is None:
                    logging.error(
                        "Could not infer atoms per frame for velocities; skipping velocity extraction."
                    )
                    vel_pbar.close()
                    return
                else:
                    logging.info(f"Inferred atoms_per_frame={inferred}")
                    atoms = inferred
            else:
                atoms = num_atoms_per_frame

            try:
                for idx, frame in stream_velocity_frames(
                    args.velocities,
                    atoms_per_frame=atoms,
                    chunk_size=args.chunk_size_mb * 1024 * 1024,
                ):
                    summary["vel_scanned"] += 1
                    vel_pbar.update(1)
                    if idx in want_frames:
                        extraction_index = extraction_map[idx]
                        out_sub = os.path.join(args.output, f"frame_{extraction_index}")
                        os.makedirs(out_sub, exist_ok=True)
                        file_path = os.path.join(
                            out_sub, f"frame_{extraction_index}.vel"
                        )
                        with open(file_path, "w") as f:
                            f.write(f"{len(frame)}\n")
                            f.write(f"Frame {idx}\n")
                            for v in frame:
                                f.write(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n")
                        summary["vel_extracted"] += 1
            except Exception as e:
                logging.exception(f"Error while streaming velocities: {e}")
                summary["errors"].append(str(e))
            finally:
                vel_pbar.close()

        # Start workers
        coords_thread = Thread(target=coords_worker, daemon=True)
        vel_thread = Thread(target=vel_worker, daemon=True)

        coords_thread.start()
        vel_thread.start()

        # Wait for completion
        coords_thread.join()
        vel_thread.join()

    except Exception as e:
        logging.exception(f"Fatal error: {e}")
        summary["errors"].append(str(e))

    finally:
        elapsed = time.time() - start_time
        # Print/write summary
        logging.info("--- Processing Summary ---")
        logging.info(f"Coords frames scanned: {summary['coords_scanned']}")
        logging.info(f"Coords frames extracted: {summary['coords_extracted']}")
        logging.info(f"Velocity frames scanned: {summary['vel_scanned']}")
        logging.info(f"Velocity frames extracted: {summary['vel_extracted']}")
        logging.info(f"Elapsed time: {elapsed:.2f} seconds")
        if summary["errors"]:
            logging.error(f"Errors encountered ({len(summary['errors'])}):")
            for err in summary["errors"]:
                logging.error(err)
        else:
            logging.info("No errors encountered.")

        print("Processing complete. Frames have been saved to the output directory.")
