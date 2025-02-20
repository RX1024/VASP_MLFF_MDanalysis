import re
import numpy as np
import matplotlib.pyplot as plt

# 打开OUTCAR文件
outcar_file = "OUTCAR"  # 修改为OUTCAR文件的实际路径

# 提取矩阵的第一列数据
matrix_data = []
with open(outcar_file, 'r') as f:
    for line in f:
       	match = re.search(r'enthalpy is ML TOTEN\s*=\s*([\d.-]+)', line)
       	if match:
            etotal = float(match.group(1))
            matrix_data.append([etotal, 0])  # 第二列的初始值设为0

# 设置矩阵的第二列
for i, row in enumerate(matrix_data):
    row[1] = i  # 公差为1的等差数列

# 转换为NumPy数组
matrix_array = np.array(matrix_data)

# 提取数据
x_data = matrix_array[:, 1]
y_data = matrix_array[:, 0]

plt.figure(figsize=(10, 6))
plt.scatter(x_data, y_data, s = 1)
plt.legend(['Total energy'])
plt.xlabel('Ionic Step')
plt.ylabel('Energy / eV')


# 调整布局，绘制图形并保存为jpg
plt.tight_layout()
plt.savefig('energy_graph.jpg')
plt.show()
