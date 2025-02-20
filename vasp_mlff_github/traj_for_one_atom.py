import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def read_xdatcar(atom_index):
    coordinates = []
    a = 0
    with open("XDATCAR", "r") as file:
        lines = file.readlines()
        for i in range(len(lines)):
            line = lines[i].strip()
            if "IL_liquid" in line:
                a = float(lines[i + 2].split()[0])
            elif "Direct configuration" in line:
                if a == 0:
                    raise ValueError("Cell parameter 'a' not found before coordinates.")
                coord_line = lines[i + atom_index].split()
                coordinates.append([float(coord) * a for coord in coord_line[:3]])
    return coordinates

def write_to_file(coordinates):
    with open("traj.txt", "w") as file:
        for coord in coordinates:
            file.write(f"{coord[0]} {coord[1]} {coord[2]}\n")

def plot_trajectory(coordinates):
    x, y, z = zip(*coordinates)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # 创建一个颜色渐变列表
    colors = plt.cm.jet(np.linspace(0, 1, len(coordinates)))

    # 绘制点和线
    for i in range(len(coordinates) - 1):
        ax.plot([x[i], x[i+1]], [y[i], y[i+1]], [z[i], z[i+1]], color=colors[i])
        ax.scatter(x[i], y[i], z[i], color=colors[i])

    # 设置坐标轴范围
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 24)
    ax.set_zlim(0, 24)

    # 设置坐标轴标签
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')

    # 添加颜色条
    mappable = plt.cm.ScalarMappable(cmap=plt.cm.jet)
    mappable.set_array(range(len(coordinates)))
    plt.colorbar(mappable, ax=ax, orientation='vertical', label='Trajectory Order')

    plt.title('3D Trajectory of Atom with Color Gradient')
    plt.show()

def main():
    atom_index = int(input("Enter the atom index: "))
    coordinates = read_xdatcar(atom_index)
    write_to_file(coordinates)
    plot_trajectory(coordinates)
    print("Coordinates written to traj.txt and plot generated.")

if __name__ == "__main__":
    main()
