import numpy as np

def read_lattice_vectors_and_coordinates(file_path, skip_frames):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    name = lines[0].strip()
    lattice_vectors = []
    coordinates = []
    current_lattice = []
    frame_counter = -1  # Initialize to -1 to account for the first frame being read
    for line in lines[2:]:  # Skip the first two lines
        if line.strip() == name:
            frame_counter += 1
            if frame_counter % skip_frames != 0:
                continue
            if current_lattice:
                lattice_vectors.append(current_lattice)
            current_lattice = []
            coordinates.append([])
        elif len(current_lattice) < 3 and len(line.split()) == 3:
            current_lattice.append([float(x) for x in line.split()])
        elif 'Direct configuration' in line:
            continue
        elif coordinates and len(line.split()) == 3:
            coordinates[-1].append([float(x) for x in line.split()])

    if current_lattice:
        lattice_vectors.append(current_lattice)

    return np.array(lattice_vectors), np.array(coordinates)

def convert_to_real_coordinates(lattice_vectors, coordinates):
    real_coordinates = []
    for i, frame in enumerate(coordinates):
        real_frame = np.dot(frame, lattice_vectors[i % len(lattice_vectors)])
        real_coordinates.append(real_frame)
    return np.array(real_coordinates)

def calculate_velocity(real_coordinates, time_step, atom_indices):
    velocities = (real_coordinates[1:, atom_indices, :] - real_coordinates[:-1, atom_indices, :]) / time_step
    return velocities

def calculate_vacf(velocities):
    n_frames, n_atoms, _ = velocities.shape
    vacf = np.zeros(n_frames-1)
    for i in range(n_atoms):  # 遍历所有原子
        for t in range(1, n_frames):  # 对于每个时间延迟t
            vel_correlations = [np.dot(velocities[tau, i], velocities[tau + t, i]) 
                                for tau in range(n_frames - t)]  # 在所有可能的时间tau上进行速度自相关
            vacf[t-1] += np.mean(vel_correlations)  # 平均所有时间点上的自相关值
    vacf /= n_atoms  # 对所有原子求平均
    if vacf[0] != 0:
        vacf /= vacf[0]
    return vacf

def parse_atom_indices(atom_indices_str):
    indices = []
    for part in atom_indices_str.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            indices.extend(range(start - 1, end))
        else:
            indices.append(int(part) - 1)
    return indices

# 用户输入
atom_indices_str = input("Enter the atom indices (e.g., '1-5' or '1,2,3' or '1-5,6,7'): ")
skip_frames = int(input("Enter the number of frames to skip between reads: "))  # Now it is the number of frames to skip
real_time_step = float(input("Enter the real time step between frames in femtoseconds: "))  # Real time step in femtoseconds
file_path = input("Enter the file path for XDATCAR: ")


# 处理读取的数据
lattice_vectors, fractional_coordinates = read_lattice_vectors_and_coordinates(file_path, skip_frames)
# 将坐标从埃转换为米
real_coordinates = convert_to_real_coordinates(lattice_vectors, fractional_coordinates) * 1e-10

atom_indices = parse_atom_indices(atom_indices_str)
# 将时间步长从飞秒转换为秒
real_velocities = calculate_velocity(real_coordinates, real_time_step * skip_frames * 1e-15, atom_indices)
vacf = calculate_vacf(real_velocities)

# 将 VACF 数据保存到文件
time_array = np.arange(len(vacf)) * real_time_step * skip_frames  # 时间以实际间隔（考虑跳帧）计算，仍然以飞秒为单位
np.savetxt("vacf_data.txt", np.column_stack((time_array, vacf)), header="Time (fs) VACF", fmt='%e')

# 假设对 VACF 进行积分以估算扩散系数
# 由于速度是以m/s为单位的，所以积分结果直接就是m^2/s，不需要再转换
D = np.trapz(vacf, dx=real_time_step * skip_frames * 1e-15) / 3

print("Self-Averaged Velocity Auto-Correlation Function (VACF) and Time have been saved to 'vacf_data.txt'")
print("Estimated Diffusion Coefficient (D):", D, "m^2/s")
