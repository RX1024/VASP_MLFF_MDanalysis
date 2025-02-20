import matplotlib.pyplot as plt
import numpy as np

# =============== 解析原子索引字符串 ===============
def parse_atom_indices(input_str):
    """
    将类似 "5,13,21" 或 "25-30" 的字符串解析成原子索引列表。
    """
    atoms = set()
    for part in input_str.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            atoms.update(range(start, end + 1))
        else:
            atoms.add(int(part))
    return sorted(atoms)

# =============== 最小镜像条件处理周期性边界 ===============
def apply_minimum_image(coord1, coord2, box_length):
    """
    对坐标差值应用最小镜像条件。
    """
    delta = np.array(coord1) - np.array(coord2)
    delta -= np.round(delta / box_length) * box_length
    return delta

# =============== 读取 XDATCAR 文件并提取坐标和晶胞大小 ===============
def read_xdatcar(oxygen_indices, hydrogen_indices):
    """
    读取 XDATCAR 文件，返回:
    - coordinates: dict, key 为原子索引，value 为每一帧的坐标 list
    - box_lengths: list, 每帧晶胞大小
    """
    coordinates = {atom_index: [] for atom_index in oxygen_indices + hydrogen_indices}
    box_lengths = []

    with open("XDATCAR", "r") as file:
        lines = file.readlines()
        for i in range(len(lines)):
            line = lines[i].strip()
            if "Direct configuration" in line:
                # 提取晶胞边长（向前找 5 行）
                scaling_factor = float(lines[i - 6].strip())
                cell_x = float(lines[i - 5].split()[0]) * scaling_factor
                cell_y = float(lines[i - 4].split()[1]) * scaling_factor
                cell_z = float(lines[i - 3].split()[2]) * scaling_factor
                box_length = (cell_x + cell_y + cell_z) / 3.0  # 假设晶胞为立方体
                box_lengths.append(box_length)

                # 解析分数坐标并转换为实际坐标
                for atom_index in coordinates:
                    coord_line = lines[i + atom_index].split()
                    fractional_coords = [float(coord) for coord in coord_line[:3]]
                    actual_coords = [fc * box_length for fc in fractional_coords]
                    coordinates[atom_index].append(actual_coords)

    return coordinates, box_lengths

# =============== 查找某个氧原子所配对的两个氢原子 ===============
def find_hydrogen_for_oxygen(coordinates, oxygen_index, hydrogen_indices, step, box_length):
    """
    对给定 step（帧）下的氧原子，寻找距离 < 1.5 Å 的两个氢原子作为配对。
    """
    oxygen_coord = coordinates[oxygen_index][step]
    hydrogen_coords = [coordinates[h_index][step] for h_index in hydrogen_indices]
    hydrogen_distances = [
        np.linalg.norm(apply_minimum_image(h_coord, oxygen_coord, box_length)) 
        for h_coord in hydrogen_coords
    ]
    # 选取距离最近的两颗氢
    hydrogen_pair = sorted(zip(hydrogen_indices, hydrogen_distances), key=lambda x: x[1])[:2]
    if hydrogen_pair[0][1] < 1.5 and hydrogen_pair[1][1] < 1.5:
        return [hydrogen_pair[0][0], hydrogen_pair[1][0]]
    return None

# =============== 计算 O-H 距离和 H-O-H 角度 ===============
def compute_OH_angle_and_distances(coordinates, oxygen_index, hydrogen_pair, step, box_length):
    """
    计算 H-O-H 角度和各自的 O-H 键长。
    """
    O = np.array(coordinates[oxygen_index][step])
    H1 = np.array(coordinates[hydrogen_pair[0]][step])
    H2 = np.array(coordinates[hydrogen_pair[1]][step])

    OH1 = apply_minimum_image(H1, O, box_length)
    OH2 = apply_minimum_image(H2, O, box_length)
    distance1 = np.linalg.norm(OH1)
    distance2 = np.linalg.norm(OH2)

    # 计算角度 (H1-O-H2)
    angle = np.arccos(
        np.dot(OH1, OH2) / (np.linalg.norm(OH1) * np.linalg.norm(OH2))
    ) * 180.0 / np.pi

    return angle, [distance1, distance2]

# =============== 计算所有步骤的角度和距离 ===============
def calculate_angle_and_distance(coordinates, oxygen_indices, hydrogen_indices, box_lengths):
    """
    返回:
    - all_angles: 所有帧中所有 H-O-H 角度的列表
    - all_distances: 所有帧中所有 O-H 距离的列表
    - total_angles: 角度数量（H-O-H 角的个数）
    - total_distances: 键长数量（O-H 键的个数）
    - low_distance_frames: 出现小于 1 Å 的键长的帧索引列表
    - low_distance_pairs: 出现小于 1 Å 键长的 (氧, 氢, 距离) 元组列表
    """
    all_angles = []
    all_distances = []
    total_angles = 0
    total_distances = 0
    low_distance_frames = []
    low_distance_pairs = []

    # 假设每帧原子数量不变
    num_steps = len(next(iter(coordinates.values())))

    for step in range(num_steps):
        frame_angles = []
        frame_distances = []
        box_length = box_lengths[step]

        for oxygen_index in oxygen_indices:
            hydrogen_pair = find_hydrogen_for_oxygen(
                coordinates, oxygen_index, hydrogen_indices, step, box_length
            )
            if hydrogen_pair:
                angle, distances = compute_OH_angle_and_distances(
                    coordinates, oxygen_index, hydrogen_pair, step, box_length
                )
                frame_angles.append(angle)
                frame_distances.extend(distances)

                total_angles += 1
                total_distances += 2  # 每对 O-(H1,H2) 含两个 O-H 键

                # 检查是否有小于 1 Å 的情况
                for dist, h_index in zip(distances, hydrogen_pair):
                    if dist < 1.0:
                        low_distance_frames.append(step)
                        low_distance_pairs.append((oxygen_index, h_index, dist))

        all_angles.extend(frame_angles)
        all_distances.extend(frame_distances)

    return (all_angles,
            all_distances,
            total_angles,
            total_distances,
            low_distance_frames,
            low_distance_pairs)

# =============== 绘制并保存概率分布图 ===============
def plot_probability_distribution(data,
                                 bins=50,
                                 range_limits=None,
                                 title="Distribution",
                                 xlabel="Value",
                                 plot_filename="distribution.jpg",
                                 data_filename="distribution_data.txt"):
    """
    绘制并保存概率分布图，支持自定义 range_limits:
       - 若 range_limits 为 None，则使用 data 的 min/max
       - 若 range_limits 为 (a, b)，则在 (a, b) 之间划分 bins 个 bin
    """
    # 如果未指定 range_limits，则让 np.histogram 自动使用数据最小值和最大值
    if range_limits is None:
        counts, bin_edges = np.histogram(data, bins=bins, density=True)
    else:
        counts, bin_edges = np.histogram(data, bins=bins, range=range_limits, density=True)

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    # probability = counts * np.diff(bin_edges) # 若想使柱子面积之和 = 1
    # 这里直接使用 counts 也行，但注意 counts 是概率密度，面积才是 1
    probabilities = counts * np.diff(bin_edges)

    # 绘制分布曲线
    plt.figure(figsize=(8, 6))
    plt.plot(bin_centers, probabilities, label=title)
    plt.xlabel(xlabel)
    plt.ylabel('Probability')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(plot_filename)
    plt.close()

    # 将中心点和概率保存到文件
    np.savetxt(data_filename,
               np.vstack((bin_centers, probabilities)).T,
               fmt='%10.5f',
               header='Bin_Center Probability')

# =============== 主程序 ===============
def main():
    # 1. 指定要分析的氧原子与氢原子
    oxygen_input = "4,12,20"   # 例如 O 原子索引：5,13,21
    hydrogen_input = "31-36"   # 例如 H 原子索引：25-30
    oxygen_indices = parse_atom_indices(oxygen_input)
    hydrogen_indices = parse_atom_indices(hydrogen_input)

    # 2. 读取坐标与晶胞大小
    coordinates, box_lengths = read_xdatcar(oxygen_indices, hydrogen_indices)

    # 3. 计算所有帧的角度、距离等信息
    (angle_values,
     distance_values,
     total_angles,
     total_distances,
     low_distance_frames,
     low_distance_pairs) = calculate_angle_and_distance(
        coordinates, oxygen_indices, hydrogen_indices, box_lengths
    )

    # ============ 下面这两行可根据先验知识来改动，若不想固定范围可设 None ============
    # 例如 H-O-H 角常见范围在 80~140° 之间，O-H 键长常见范围在 0.5~2.0 Å 之间
    angle_range = (112, 136)
    distance_range = (1.18, 1.42)
    # ===============================================================================

    # 4. 绘制并保存概率分布图（使用统一上下限）
    plot_probability_distribution(angle_values,
                                 bins=50,
                                 range_limits=angle_range,
                                 title='H-O-H Angle Distribution',
                                 xlabel='Angle (degrees)',
                                 plot_filename='angle_distribution.jpg',
                                 data_filename='angle_distribution_data.txt')

    plot_probability_distribution(distance_values,
                                 bins=50,
                                 range_limits=distance_range,
                                 title='O-H Distance Distribution',
                                 xlabel='Distance (Å)',
                                 plot_filename='distance_distribution.jpg',
                                 data_filename='distance_distribution_data.txt')

    # 5. 打印统计信息
    avg_angle = np.mean(angle_values)
    avg_distance = np.mean(distance_values)
    print(f"Average H-O-H angle: {avg_angle:.3f} degrees")
    print(f"Average O-H distance: {avg_distance:.3f} Å")
    print(f"Total number of angles (H-O-H): {total_angles}")
    print(f"Total number of distances (O-H): {total_distances}")


    print("Probability distribution plots generated and saved.")

# 入口
if __name__ == "__main__":
    main()
