#!/usr/bin/env python3
import os
import numpy as np

def indent_lines(text, indent=" "):
    """
    text を行分割して、各行の先頭をチェックし、
    - 行頭が '-' (マイナス) ならインデントを入れない
    - それ以外の行なら先頭に半角スペースを入れる
    """
    lines = text.splitlines()
    out_lines = []
    for line in lines:
        if line.startswith('-'):
            # 行がマイナス記号で始まる場合はインデントを入れない
            out_lines.append(line)
        else:
            # それ以外の場合はインデントを付与
            out_lines.append(indent + line)
    return "\n".join(out_lines)

def extract_hessian_custom_order(lines, atom_num):
    hessian_start = lines.index("$hessian\n") + 3
    hessian_end = hessian_start + 1 + ((atom_num*3 // 5 + 1) * atom_num*3)
    hessian_lines = lines[hessian_start : hessian_end]

    row_matrix = []
    for line in hessian_lines:
        parts = line.split()
        numbers = [float(p) for p in parts if ("E" in p or "." in p)]
        if numbers:
            row_matrix.append(numbers)

    group_size = atom_num**2
    num_groups = len(row_matrix) // group_size
    groups = [row_matrix[i * group_size : (i + 1) * group_size] for i in range(num_groups)]

    matrix = []
    for i in range(group_size):
        combined_row = []
        for g in groups:
            combined_row.extend(g[i])
        matrix.append(combined_row)

    custom_order = []
    for i, row in enumerate(matrix):
        for j in range(i + 1):
            custom_order.append(str(row[j]))

    return custom_order

def process_directory(work_dir: str, atom_num: int, line_flag: bool):
    hess_path   = os.path.join(work_dir, "job.hess")
    out_path    = os.path.join(work_dir, "job.out")
    engrad_path = os.path.join(work_dir, "job.engrad")

    os.makedirs("../minfo.files", exist_ok=True)
    minfo_path  = os.path.join("../minfo.files", f"{work_dir}.minfo")

    with open(hess_path, "r", encoding="UTF-8") as f:
        hessdata = f.readlines()

    # Dipole Derivatives
    dipole_index = hessdata.index("$dipole_derivatives\n") + 2
    dipole_line = hessdata[dipole_index : dipole_index + atom_num * 3]
    dipole_line = [line.strip().split() for line in dipole_line]
    dipole_data = []
    for col in range(3):
        for row in range(atom_num * 3):
            dipole_data.append(dipole_line[row][col])
    dipole = "\n".join(
        ", ".join(dipole_data[i : i+5])
        for i in range(0, len(dipole_data), 5)
    )

    # Hessian (custom order)
    upper_triangle_line = extract_hessian_custom_order(hessdata, atom_num)
    upper_triangle = "\n".join(
        ", ".join(upper_triangle_line[i : i+5])
        for i in range(0, len(upper_triangle_line), 5)
    )

    # out ファイルからエネルギーなどを取得
    with open(out_path, "r", encoding="UTF-8") as f:
        out_data = f.readlines()

    energy_index = next(
        i for i, line in enumerate(out_data)
        if "TOTAL SCF ENERGY" in line
    ) + 3
    energy_str = out_data[energy_index].split("Total Energy       :")[1].split("Eh")[0].strip()
    energy = float(energy_str)

    charge_index = next(i for i, line in enumerate(out_data) if "Total Charge" in line)
    charge_str = out_data[charge_index].split("....")[1].strip()
    mult_str   = out_data[charge_index + 1].split("....")[1].strip()
    charge = float(charge_str)
    mult   = float(mult_str)

    di_mo_index = next(i for i, line in enumerate(out_data) if "Total Dipole Moment" in line)
    dipole_moment = out_data[di_mo_index].split(":")[1].strip().split()
    di_moment = ", ".join(dipole_moment)

    pol_index = next(i for i, line in enumerate(out_data) if "The raw cartesian tensor (atomic units):" in line)
    pol_row = out_data[pol_index + 1 : pol_index + 4]
    pol_line = list(range(6))
    for i, line_ in enumerate(pol_row):
        parts = line_.split()
        if i == 0:
            pol_line[0] = parts[0]
            pol_line[1] = parts[1]
            pol_line[2] = parts[2]
        elif i == 1:
            pol_line[3] = parts[1]
            pol_line[4] = parts[2]
        elif i == 2:
            pol_line[5] = parts[2]
    pol = "\n".join(
        ", ".join(str(x) for x in pol_line[i : i+5])
        for i in range(0, len(pol_line), 5)
    )

    with open(engrad_path, "r", encoding="UTF-8") as f:
        g_data = f.readlines()
    grad_index = next(i for i, line in enumerate(g_data) if "The current gradient" in line)
    grad_line = g_data[grad_index + 2 : grad_index + 2 + atom_num * 3]
    grad_line = [line.strip() for line in grad_line]
    grad = "\n".join(
        ", ".join(grad_line[i : i+5])
        for i in range(0, len(grad_line), 5)
    )

    with open(minfo_path, "w", encoding="UTF-8") as f:
        f.write("# minfo File version 2:\n")
        f.write("#\n")
        f.write("[ Electronic Data ]\n")

        f.write("Energy\n")
        f.write(f"{energy}\n")

        f.write("Charge\n")
        f.write(f"{charge}\n")

        f.write("Multiplicity\n")
        f.write(f"{mult}\n")

        f.write("Gradient\n")
        grad_count_str = str(atom_num*3)
        f.write(grad_count_str + "\n")
        # 勾配の各行の先頭文字が '-' なら空白なし、その他は空白を入れる
        f.write(indent_lines(grad) + "\n")

        f.write("Hessian\n")
        hessian_count = atom_num*3 * (atom_num*3 + 1) // 2
        f.write(str(hessian_count) + "\n")
        f.write(indent_lines(upper_triangle) + "\n")

        f.write("Dipole Moment\n")
        f.write("3\n")
        f.write(indent_lines(di_moment) + "\n")

        f.write("Polarizability\n")
        f.write("6\n")
        f.write(indent_lines(pol) + "\n")

        f.write("Dipole Derivative\n")
        dipole_count = atom_num*3*3
        f.write(str(dipole_count) + "\n")
        f.write(indent_lines(dipole.replace("E","e")) + "\n")
        f.write("\n")
    print(f"[DONE] {work_dir}:{work_dir}.minfo を出力しました.")

def main():
    atom_num_str = input("原子数を入力してください（例：3）: ")
    line_flag_str = input("直線分子かどうか True or False で入力してください: ")

    try:
        atom_num = int(atom_num_str)
    except:
        print("原子数は整数で入力してください。例: 3")
        return
    
    if line_flag_str.lower() == "true":
        line_flag = True
    elif line_flag_str.lower() == "false":
        line_flag = False
    else:
        print("直線分子フラグは True か False を入力してください。")
        return

    current_dir = os.getcwd()
    for d in sorted(os.listdir(current_dir)):
        if os.path.isdir(d) and d.startswith("mkqff"):
            print(f"[INFO] 処理ディレクトリ: {d}")
            process_directory(d, atom_num, line_flag)

if __name__ == "__main__":
    main()

