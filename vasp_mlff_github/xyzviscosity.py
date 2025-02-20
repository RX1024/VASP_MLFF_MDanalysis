import numpy as np
import matplotlib.pyplot as plt

def parse_stress_and_volume_from_outcar(filename="OUTCAR"):
    """从OUTCAR文件中解析应力和体积数据"""
    stress_data = []
    volume_data = []
    recent_volume = None
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "volume of cell :" in line:
                recent_volume = float(line.split()[-1])
            if "in kB" in line and recent_volume is not None:
                stress = line.split()[2:8]
                stress = [float(s) for s in stress]
                stress_data.append(stress)
                volume_data.append(recent_volume)
                recent_volume = None
    return np.array(stress_data), np.array(volume_data)

def calculate_acf(stress_component):
    """计算自相关函数"""
    num_points = len(stress_component)
    acf = np.correlate(stress_component, stress_component, mode='full') / num_points
    acf = acf[num_points-1:]  # 只取正的时间延迟部分
    return acf

def integrate_acf(acf, dt):
    """数值积分自相关函数以估算粘度"""
    return np.trapz(acf, dx=dt)

def calculate_viscosity_over_time(stress_data, volume_data, temperature, dt):
    """计算并返回每个时间步的粘度"""
    kB = 1.380649e-23  # 玻尔兹曼常数，单位J/K
    dt = dt * 1e-15  # 时间步长从飞秒转换为秒
    viscosities_over_time = []

    for i, direction in enumerate(['XY', 'YZ', 'ZX']):
        stress_component = stress_data[:, i+3] * 1e8  # 千巴转换为Pascal
        acf = calculate_acf(stress_component)
        viscosity_contributions = []
        for j in range(len(acf)):
            viscosity = integrate_acf(acf[:j+1], dt) * volume_data[j] * 1e-30 / (kB * temperature)
            viscosity_contributions.append(viscosity)
        viscosities_over_time.append(viscosity_contributions)

    # 计算平均粘度
    average_viscosity = np.mean(viscosities_over_time, axis=0)
    viscosities_over_time.append(average_viscosity)
    return viscosities_over_time

def plot_viscosity_over_time(viscosities_over_time, dt):
    """绘制粘度随时间变化的图"""
    time_steps = np.arange(len(viscosities_over_time[0])) * dt * 1e-15
    plt.figure(figsize=(10, 6))
    directions = ['XY', 'YZ', 'ZX', 'Average']
    for i, direction in enumerate(directions):
        plt.plot(time_steps, viscosities_over_time[i], label=f'Viscosity in {direction} direction')
    plt.xlabel('Time (s)')
    plt.ylabel('Viscosity (Pa·s)')
    plt.title('Viscosity Over Time')
    plt.legend()
    plt.show()

def save_viscosity_data(viscosities_over_time, filename="viscosities.txt"):
    """保存粘度数据到文本文件"""
    np.savetxt(filename, np.transpose(viscosities_over_time), fmt='%.5e',
               header='Viscosity in XY, Viscosity in YZ, Viscosity in ZX, Average Viscosity')

if __name__ == "__main__":
    outcar_file = "OUTCAR"
    temperature = 300  # 单位: K
    dt = 1000  # MD时间步长，单位: fs

    stress_data, volume_data = parse_stress_and_volume_from_outcar(outcar_file)
    viscosities_over_time = calculate_viscosity_over_time(stress_data, volume_data, temperature, dt)
    plot_viscosity_over_time(viscosities_over_time, dt)
    save_viscosity_data(viscosities_over_time)
