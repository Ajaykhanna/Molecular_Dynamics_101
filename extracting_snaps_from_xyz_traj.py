import os
import sys
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt


# Filepath for the large .xyz file
file_path = str(sys.argv[1])  # 'coords.xyz'

def extract_frames_as_xyz(trajectory_data, start=0, stop=None, step=1, output_dir="extracted_frames"):
    """
    Extracts specific frames from the trajectory and saves them as individual .xyz files.
    
    Parameters:
        trajectory_data (list): List of frames where each frame is a list of atom data.
        start (int): Starting frame index (default is 0).
        stop (int): Stopping frame index (default is len(trajectory_data)).
        step (int): Step size for extraction (default is 1).
        output_dir (str): Directory to save the extracted frames (default is 'extracted_frames').
    """
    # Set default stop value if not provided
    if stop is None:
        stop = len(trajectory_data)
    
    # Validate inputs
    if start < 0 or stop > len(trajectory_data) or step <= 0:
        raise ValueError("Invalid start, stop, or step values for frame extraction.")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract frames based on user-specified range
    frame_indices = range(start, stop, step)
    for idx, frame_num in enumerate(frame_indices):
        frame = trajectory_data[frame_num]
        file_path = os.path.join(output_dir, f"frame_{idx + 1}.xyz")
        with open(file_path, "w") as f:
            f.write(f"{len(frame)}\n")
            f.write(f"Frame {frame_num}\n")
            for atom in frame:
                f.write(f"{atom[0]} {atom[1]:.10f} {atom[2]:.10f} {atom[3]:.10f}\n")

    return output_dir


# Optimized version of the dihedral angle calculation function
def find_dihedral_angle(a, b, c, d):
    ba = a - b
    bc = c - b
    cd = d - c
    b_norm = ba / np.linalg.norm(ba)
    c_norm = bc / np.linalg.norm(bc)
    d_norm = cd / np.linalg.norm(cd)
    x = np.cross(b_norm, c_norm)
    y = np.cross(c_norm, d_norm)
    cos_angle = np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y))
    cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Ensure no domain errors
    angle = np.degrees(np.arccos(cos_angle))
    return angle - 180.0

# Integrating this functionality into a dihedral analysis script
def calculate_dihedral_angles(extracted_data):
    dihedral_angles = []
    
    for i, frame in enumerate(extracted_data):
        # Extract relevant atoms for dihedral calculation
        try:
            a = np.array([frame[4][1], frame[4][2], frame[4][3]])
            b = np.array([frame[1][1], frame[1][2], frame[1][3]])
            c = np.array([frame[0][1], frame[0][2], frame[0][3]])
            d = np.array([frame[2][1], frame[2][2], frame[2][3]])
            dihedral_angle = find_dihedral_angle(a, b, c, d)
            dihedral_angles.append(dihedral_angle)
        except IndexError:
            print(f"Error in frame {i}: Missing required atoms for dihedral calculation.")
            continue

    # Save dihedral angles to a file
    output_file = "dihedral_angles_traj.txt"
    with open(output_file, "w") as f:
        for i, angle in enumerate(dihedral_angles):
            f.write(f"Frame {i}: {angle:.2f} degrees\n")

    # Plot the distribution of dihedral angles
    plt.hist(dihedral_angles, bins=50, edgecolor='k', alpha=0.7)
    plt.title("Distribution of Dihedral Angles")
    plt.xlabel("Dihedral Angle (degrees)")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.savefig("dihedral_angle_distribution.png")
    #plt.show()

    return output_file, "dihedral_angle_distribution.png"


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

        frame_data = lines[i:i + frame_size]
        frame_results = []
        for atom_line in frame_data[2:]:  # Skip the first two lines (header + metadata)
            parts = atom_line.split()
            if len(parts) == 4:  # Ensure line has atomic symbol and 3 coordinates
                atom_symbol, x, y, z = parts
                frame_results.append((atom_symbol, float(x), float(y), float(z)))

        results.append(frame_results)
        i += frame_size

    return results, leftover

# Function to split the file into manageable chunks for multiprocessing
def chunk_reader(file_path, chunk_size=1024 * 1024 * 50):  # 50MB per chunk
    with open(file_path, 'r') as file:
        leftover = ""
        while True:
            chunk = file.read(chunk_size)
            if not chunk:
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

# Main script
if __name__ == "__main__":
    extracted_data = process_file_in_parallel(file_path)
    print(f"Extracted frames: {len(extracted_data)}")
    dihedral_output, plot_output = calculate_dihedral_angles(extracted_data)
    output_frames_dir = extract_frames_as_xyz(extracted_data, start=100000, stop=100001, step=1)

