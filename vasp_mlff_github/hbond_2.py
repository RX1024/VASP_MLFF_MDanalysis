import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import random


def parse_range(input_str):
    ranges = input_str.split(',')
    result = []
    for r in ranges:
        if '-' in r:
            start, end = map(int, r.split('-'))
            result.extend(range(start - 1, end))  # 转换为从0开始的索引
        else:
            result.append(int(r) - 1)
    return result


def get_name_line(file_path):
    with open(file_path, 'r') as file:
        name_line = file.readline().strip()
    return name_line


def read_frame(file, name_line):
    frame = []
    line = file.readline()
    while line:
        if 'Cartesian' in line:
            break
        line = file.readline()
    line = file.readline()
    while line:
        if name_line in line:
            break
        frame.append(list(map(float, line.split())))
        line = file.readline()
    return frame


def calculate_distance(coord1, coord2):
    return np.linalg.norm(np.array(coord1) - np.array(coord2))


def calculate_angle(coordA, coordB, coordH):
    BA = np.array(coordA) - np.array(coordB)
    BH = np.array(coordH) - np.array(coordB)

    if np.isclose(np.linalg.norm(BA), 0) or np.isclose(np.linalg.norm(BH), 0):
        return np.nan  # 返回NaN，表示角度无法计算

    cosine_angle = np.dot(BA, BH) / (np.linalg.norm(BA) * np.linalg.norm(BH))
    cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)


def find_hydrogen_bonds(frame, A_indices, B_indices, H_indices, max_length, max_angle, H_A_max_distance=1.5):
    h_bonds = []
    valid_H_for_A = {A: [] for A in A_indices}
    
    for H in H_indices:
        for A in A_indices:
            if H != A:
                dist_HA = calculate_distance(frame[H], frame[A])
                if dist_HA < H_A_max_distance:
                    valid_H_for_A[A].append((H, dist_HA))

    for B in B_indices:
        for A in A_indices:
            if B != A:
                for H, dist_HA in valid_H_for_A[A]:
                    if H != B:
                        dist_BH = calculate_distance(frame[B], frame[H])
                        if dist_BH < max_length:
                            angle = 180 - calculate_angle(frame[A], frame[H], frame[B])
                            if angle < max_angle:
                                h_bonds.append((dist_BH, angle, A, B, H))

    return h_bonds


def plot_hydrogen_bond_distribution(h_bond_lengths, h_bond_angles, max_length, max_angle):
    plt.rcParams['font.family'] = 'Calibri'
    plt.rcParams['font.size'] = 16
    plt.rcParams['axes.linewidth'] = 1.8
    xy = np.vstack([h_bond_lengths, h_bond_angles])
    kde = gaussian_kde(xy)
    z = kde(xy)  # 直接对数据点计算密度

    # 筛选 <3Å 和 <100° 范围内的点
    valid_indices = (np.array(h_bond_lengths) < 3) & (np.array(h_bond_angles) < 100)
    valid_point_count = np.sum(valid_indices)  # 计算范围内点的数量
    print(f"在指定范围 (<3Å, <100°) 内的点数量：{valid_point_count}")
    
    if np.any(valid_indices):  # 确保范围内有点
        valid_lengths = np.array(h_bond_lengths)[valid_indices]
        valid_angles = np.array(h_bond_angles)[valid_indices]
        valid_densities = z[valid_indices]

        # 找到范围内密度最大点
        max_idx = np.argmax(valid_densities)
        x_peak = valid_lengths[max_idx]
        y_peak = valid_angles[max_idx]
        density_peak = valid_densities[max_idx]

        print("密度最大点及密度值（长度, 角度, 密度）：")
        print(f"Length: {x_peak:.2f}, Angle: {y_peak:.2f}, Density: {density_peak:.2e}")
    else:
        print("在指定范围 (<3Å, <100°) 内未找到有效点。")
        x_peak, y_peak, density_peak = None, None, None

    # 绘制分布图
    fig, ax = plt.subplots()

    # 使用自适应颜色映射范围，自动根据密度分布选择最小和最大值
    scatter = ax.scatter(h_bond_lengths, h_bond_angles, c=z, s=20, edgecolor='none', 
                         cmap='coolwarm')
    
    # 添加颜色条（自动根据数据自适应最大值和最小值）
    cbar = plt.colorbar(scatter, ax=ax, label='Density')
    cbar.ax.tick_params(labelsize=16, width=3, length=5)
    cbar.set_label('Density', fontweight='bold')

    # 如果需要突出显示最大密度点，可以取消下面的注释
    # if x_peak is not None and y_peak is not None:
    #     ax.plot(x_peak, y_peak, 'go')  # 使用绿色标记最大密度点
    #     ax.text(x_peak, y_peak, f"  ({x_peak:.2f}, {y_peak:.2f})", fontsize=14, color='green', fontweight='bold')

    ax.set_xlabel('Length (Å)', fontweight='bold')
    ax.set_ylabel('Angle (°)', fontweight='bold')
    ax.set_xlim(1.5, max_length)
    ax.set_ylim(0, max_angle)
    ax.tick_params(axis='both', labelsize=16, width=3, length=5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    cbar.ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    
    # 设置轴上的数字加粗
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight('bold')

    # 保存图像为文件
    image_filename = "hydrogen_bond_distribution_optimized.jpg"
    plt.savefig(image_filename, dpi=300, bbox_inches='tight')
    print(f"氢键分布图已保存为文件：{image_filename}")
    plt.show()


def main():
    name_line = get_name_line('XDATCAR_EX')
    A = input("请输入供体原子序号数(A)，如3，4，5或1-5: ")
    B = input("请输入受体原子序号数(B)，如3，4，5或1-5: ")
    H = input("请输入氢原子序号数(H)，如3，4，5或1-5: ")
    max_length = float(input("请输入氢键最大长度（L）: "))
    max_angle = float(input("请输入氢键最大角度（D）: "))

    A_indices = parse_range(A)
    B_indices = parse_range(B)
    H_indices = parse_range(H)

    with open('XDATCAR_EX', 'r') as file:
        h_bond_lengths = []
        h_bond_angles = []
        frame_number = 0
        while True:
            frame = read_frame(file, name_line)
            if not frame:
                break
            frame_number += 1
            h_bonds = find_hydrogen_bonds(frame, A_indices, B_indices, H_indices, max_length, max_angle)
            h_bond_lengths.extend([bond[0] for bond in h_bonds])
            h_bond_angles.extend([bond[1] for bond in h_bonds])

            # 输出每帧找到的氢键数量
            print(f"Frame {frame_number}: {len(h_bonds)} hydrogen bonds found.")

    total_h_bonds = len(h_bond_lengths)
    print(f"Total hydrogen bonds found: {total_h_bonds}")

    # 输出符合条件的点的数量
    valid_indices = (np.array(h_bond_lengths) < 3) & (np.array(h_bond_angles) < 100)
    valid_point_count = np.sum(valid_indices)
    print(f"符合条件的氢键数量 (<3Å, <100°)：{valid_point_count}")

    r = int(input("请输入随机选取的氢键数目 (r): "))
    if r > total_h_bonds:
        print(f"输入值超过总氢键数，将使用所有 {total_h_bonds} 个氢键。")
        r = total_h_bonds

    selected_indices = random.sample(range(total_h_bonds), r)
    selected_h_bond_lengths = [h_bond_lengths[i] for i in selected_indices]
    selected_h_bond_angles = [h_bond_angles[i] for i in selected_indices]

    plot_hydrogen_bond_distribution(selected_h_bond_lengths, selected_h_bond_angles, max_length, max_angle)


if __name__ == "__main__":
    main()
