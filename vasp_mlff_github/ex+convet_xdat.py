def remove_last_line(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # 仅当文件不为空时才移除最后一行
    if lines:
        lines = lines[:-1]

    with open(file_path, 'w') as file:
        file.writelines(lines)
def get_first_line(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
    return first_line

def get_first_line(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
    return first_line

def count_name_lines(file_path, name_line):
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() == name_line:
                count += 1
    return count

def extract_frames(file_path, start_frame, end_frame, name_line):
    frame_count = 0
    extract = False
    extracted_lines = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() == name_line:
                if extract and frame_count == end_frame:
                    extracted_lines.append(line)  # 确保包含结束帧的name行
                    break
                frame_count += 1
                if frame_count == start_frame:
                    extract = True

            if extract:
                extracted_lines.append(line)

    with open('XDATCAR_EX', 'w') as new_file:
        new_file.writelines(extracted_lines)

def process_frames(file_path, name_line):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    frame_start = False
    for i, line in enumerate(lines):
        if line.strip() == name_line:
            frame_start = True
            # 检查索引是否在范围内
            if i + 2 < len(lines):
                side = float(lines[i + 2].split()[0])  # 假设第一个数字是float
            else:
                break  # 跳出循环，因为无法处理不完整的帧
            continue

        if frame_start and 'Direct' in line:
            lines[i] = line.replace('Direct', 'Cartesian')
            for j in range(i + 1, len(lines)):
                if j >= len(lines) or lines[j].strip() == name_line or 'Direct' in lines[j]:
                    break
                coords = lines[j].split()
                coords = [float(coord) * side for coord in coords]
                lines[j] = ' '.join(map(str, coords)) + '\n'

    with open(file_path, 'w') as file:
        file.write(name_line + '\n')  # 写回原始的第一行
        file.writelines(lines)

    remove_first_line_from_XDATCAR_EX()

def remove_first_line_from_XDATCAR_EX():
    with open('XDATCAR_EX', 'r') as file:
        lines = file.readlines()
    
    with open('XDATCAR_EX', 'w') as file:
        file.writelines(lines[1:])
# 使用示例
def main():
    file_path = 'XDATCAR'
    name_line = get_first_line(file_path)
    name_count = count_name_lines(file_path, name_line)
    print(f"文件中有 {name_count} 个 '{name_line}' 行")

    # 提取特定帧
    start_frame = int(input("请输入开始帧数: "))
    end_frame = int(input("请输入结束帧数: "))
    extract_frames(file_path, start_frame, end_frame, name_line)

    # 处理XDATCAR_EX文件
    process_frames('XDATCAR_EX', name_line)

    # 删除XDATCAR_EX文件的最后一行
    remove_last_line('XDATCAR_EX')


if __name__ == "__main__":
    main()