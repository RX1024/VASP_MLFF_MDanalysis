def extract_frames(xdatcar_filename, start_frame, end_frame, step, output_filename):
    """
    提取XDATCAR中特定范围和特定间隔的帧及其内容。

    :param xdatcar_filename: XDATCAR 文件的路径。
    :param start_frame: 起始帧（从1开始计数）。
    :param end_frame: 结束帧（不包括此帧）。
    :param step: 提取帧的间隔。
    :param output_filename: 输出文件的名称。
    """
    with open(xdatcar_filename, 'r') as file:
        lines = file.readlines()

    start_indices = [i for i, line in enumerate(lines) if "IL_liquid" in line]
    if end_frame > len(start_indices):
        raise ValueError("指定的结束帧超出了文件中的帧数。")

    # 为了提取完整的帧，需要找到每个帧的结束索引
    end_indices = start_indices[1:] + [len(lines)]

    # 根据步长提取帧
    selected_frames = range(start_frame - 1, min(end_frame, len(start_indices)), step)

    # 写入新文件
    with open(output_filename, 'w') as output_file:
        for frame in selected_frames:
            for line in lines[start_indices[frame]:end_indices[frame]]:
                output_file.write(line)

    print(f"已从 {xdatcar_filename} 中提取第 {start_frame} 到第 {end_frame-1} 帧（步长为 {step}），并保存到 {output_filename}")

# 用户输入的起始帧、结束帧和步长
start_frame = int(input("请输入起始帧（从1开始计数）: "))
end_frame = int(input("请输入结束帧（不包括此帧）: "))
step = int(input("请输入提取帧的间隔步长: "))

# XDATCAR 文件路径和输出文件名
xdatcar_filename = 'XDATCAR'  # 根据实际路径修改
output_filename = 'extracted_frames_XDATCAR'  # 可以根据需要修改

# 调用函数提取帧
extract_frames(xdatcar_filename, start_frame, end_frame, step, output_filename)
