import numpy as np
from scipy.ndimage import gaussian_filter1d

def read_lattice_vectors_and_coordinates(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    name = lines[0].strip()  # 读取第一行作为name
    lattice_vectors = []
    coordinates = []
    current_lattice = []
    for line in lines[0:]:  # 从第二行开始处理
        if line.strip() == name:  # 新一帧的开始
            if current_lattice:  # 如果之前有晶格系数，添加到列表
                lattice_vectors.append(current_lattice)
            current_lattice = []
            coordinates.append([])  # 开始新一帧的坐标收集
        elif len(current_lattice) < 3 and len(line.split()) == 3:
            current_lattice.append([float(x) for x in line.split()])
        elif 'Direct configuration' in line and coordinates:
            continue  # 找到坐标开始的标志
        elif coordinates[-1] is not None and len(line.split()) == 3:
            coordinates[-1].append([float(x) for x in line.split()])

    if current_lattice:  # 添加最后一帧的晶格系数
        lattice_vectors.append(current_lattice)

    return np.array(lattice_vectors), np.array(coordinates)

def calculate_distances(atoms, atom_indices, cell):
    distances = []

    for i in atom_indices[0]:
        for j in atom_indices[1]:
            if i != j:  # 排除自身
                diff = atoms[i] - atoms[j]
                diff = diff - np.round(diff)  # Apply PBC
                dist = np.sqrt(np.sum(np.dot(diff, cell) ** 2))  # Convert to Cartesian
                distances.append(dist)

    return np.array(distances)

def calculate_rdf(distances, num_particles, system_volume, dr=0.2, max_r=10.8):
    num_bins = int(max_r / dr)
    rdf = np.zeros(num_bins)
    bin_edges = np.linspace(0, max_r, num_bins + 1)

    for dist in distances:
        if dist > 0:  # Exclude distance to self
            bin_index = int(dist // dr)
            if bin_index < num_bins:
                rdf[bin_index] += 1

    # Normalize RDF
    density = num_particles / system_volume  # Calculate the density
    normalization = (4/3 * np.pi * (bin_edges[1:]**3 - bin_edges[:-1]**3))
    rdf /= (density * normalization)

    return rdf, bin_edges[:-1]

def parse_atom_indices(indices_str):
    if '-' in indices_str:
        start, end = map(int, indices_str.split('-'))
        return list(range(start - 1, end))
    else:
        return [int(index) - 1 for index in indices_str.split(',')]

# Path to your XDATCAR file
file_path = 'XDATCAR'

# Read lattice vectors and coordinates from file
lattice_vectors, coordinates = read_lattice_vectors_and_coordinates(file_path)

# 输入的原子序号范围
atom_indices_str1 = "193-368"  # Replace with your range
atom_indices_str2 = "193-368"  # Replace with your range
#阴：25,26,27,28,29,30,31,32,89,90,91,92,93,94,95,96
#阳：129,130,131,132,133,134,135,136,145,146,147,148,149,150,151,152
#水：193-368
# 解析原子序号
atom_indices1 = parse_atom_indices(atom_indices_str1)
atom_indices2 = parse_atom_indices(atom_indices_str2)
atom_indices = (atom_indices1, atom_indices2)

# Calculate RDF for each frame and average
all_rdfs = []
frame_count = 0  # 用于记录帧数
for frame_index, atoms in enumerate(coordinates):
    cell = np.array(lattice_vectors[frame_index])  # 使用当前帧的晶格矢量
    system_volume = np.linalg.det(cell)  # 计算当前帧的体系体积
    num_particles = len(atom_indices[0]) * len(atom_indices[1])  # 总粒子数
    distances = calculate_distances(atoms, atom_indices, cell)
    rdf, r_values = calculate_rdf(distances, num_particles, system_volume)
    all_rdfs.append(rdf)
    frame_count += 1  # 增加帧数计数

average_rdf = np.mean(all_rdfs, axis=0)

sigma = 1.5  # 高斯核的标准差，调整以改变平滑程度
smoothed_rdf = gaussian_filter1d(average_rdf, sigma)

# 打印并保存平滑后的 RDF 数据
with open('smoothed_rdf.txt', 'w') as file:
    for r, g in zip(r_values, smoothed_rdf):
        print(f'{r:.2f} {g:.3f}')
        file.write(f'{r:.2f}\t{g:.3f}\n')
# 输出处理的总帧数
print(f'Total frames processed: {frame_count}')
