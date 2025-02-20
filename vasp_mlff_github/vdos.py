import numpy as np
from numpy.fft import fft, fftfreq
import matplotlib.pyplot as plt

# 读取VACF数据
def read_vacf(file_name):
    data = np.loadtxt(file_name)
    times = data[:, 0] * 1e-15  # 时间单位：将飞秒转换为秒
    vacf_values = data[:, 1]
    return times, vacf_values

# 执行零扩展和窗函数处理
def zero_padding(vacf_values, padding_factor=1):
    padded_length = padding_factor * len(vacf_values)
    return np.pad(vacf_values, (0, padded_length - len(vacf_values)), 'constant')

# 计算VDOS
# 计算VDOS
def calculate_vdos(vacf_values, time_step):
    # 减去平均值
    vacf_values = vacf_values - np.mean(vacf_values)
    zero_padded_vacf = zero_padding(vacf_values)  # 零扩展
    vdos_freq = fft(zero_padded_vacf)  # FFT
    vdos_magnitude = np.abs(vdos_freq)  # 取模，以保证VDOS非负
    freqs = fftfreq(len(zero_padded_vacf), d=time_step)  # 频率轴
    return freqs, vdos_magnitude


# 主程序
if __name__ == '__main__':
    times, vacf_values = read_vacf("vacf_data.txt")
    time_step = times[1] - times[0]

    frequencies, vdos = calculate_vdos(vacf_values, time_step)

    # 转换频率到波数
    c = 29979245800  # 光速，单位：厘米/秒
    wave_numbers = frequencies / c
    positive_wave_numbers = wave_numbers > 0
    wave_numbers = wave_numbers[positive_wave_numbers]
    vdos = vdos[positive_wave_numbers]

    # 保存VDOS数据到文件
    np.savetxt("vdos_data.txt", np.column_stack((wave_numbers, vdos)), header="Wave Number (cm⁻¹), VDOS", fmt='%e')

    plt.plot(wave_numbers, vdos)
    plt.xlabel("Wave Number (cm⁻¹)")
    plt.ylabel("VDOS")
    plt.title("Vibrational Density of States (VDOS)")
    plt.xlim(0, max(wave_numbers))
    plt.grid(True)
    plt.show()
