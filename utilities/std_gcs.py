import sys
import numpy as np
import matplotlib.pyplot as plt


# Compute Geometric Centeroids
def compute_geometeric_center(coordinates):
    geom_center = np.mean(coordinates, axis=0)

    return geom_center


def read_file(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    return lines


def compute_difference_between_arrays(array1, array2):
    return np.linalg.norm(array1 - array2)


reference_file = str(sys.argv[1])
nAtoms = int(sys.argv[2])

lines = read_file(reference_file)
ref_coordinates = []
for line in lines[2 : nAtoms + 1]:
    ref_coordinates.append([float(x) for x in line.split()[1:]])
ref_geom_center = compute_geometeric_center(ref_coordinates)


start_frame = int(sys.argv[3])
end_frame = int(sys.argv[4])

var_geom_center = np.zeros([end_frame - start_frame + 1, 3])

for frame in range(start_frame, end_frame + 1):
    var_file = f"{frame}/solute_solvent_5.0ang_qm.xyz"
    lines = read_file(var_file)
    var_coordinates = []
    for line in lines[2 : nAtoms + 1]:
        var_coordinates.append([float(x) for x in line.split()[1:]])

    var_geom_center[frame] = compute_geometeric_center(var_coordinates)

diff = np.zeros(end_frame - start_frame + 1)

for frame in range(start_frame, end_frame + 1):
    diff[frame] = compute_difference_between_arrays(
        ref_geom_center, var_geom_center[frame]
    )

print(diff)
