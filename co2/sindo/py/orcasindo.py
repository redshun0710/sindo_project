#!/usr/bin/env python3
import os
import numpy as np

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

def process_directory(
    work_dir: str,
    atom_num: int,
    line_flag: bool,
):

    # パスを組み立て
    hess_path   = os.path.join(work_dir, f"job.hess")
    out_path    = os.path.join(work_dir, f"job.out")
    engrad_path = os.path.join(work_dir, f"job.engrad")
    minfo_path  = os.path.join(work_dir, f"job.minfo")

    with open(hess_path, "r", encoding="UTF-8") as f:
        hessdata = f.readlines()

    # 双極子微分 ($dipole_derivatives)
    dipole_index = hessdata.index("$dipole_derivatives\n") + 2
    dipole_line = hessdata[dipole_index : dipole_index + atom_num * 3]
    dipole_line = [line.strip().split() for line in dipole_line]
    dipole_data = []
    for col in range(3):
        for row in range(atom_num * 3):
            dipole_data.append(dipole_line[row][col])
    # 見やすく整形
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

    # out ファイルを読んでエネルギーなどを取得
    with open(out_path, "r", encoding="UTF-8") as f:
        out_data = f.readlines()

    energy_index = next(
        i for i, line in enumerate(out_data)
        if "*** OPTIMIZATION RUN DONE ***" in line
    ) - 3
    # 例: "FINAL SINGLE POINT ENERGY      -113.4567890"
    energy_str = out_data[energy_index].split("FINAL SINGLE POINT ENERGY")[1].strip()
    energy = float(energy_str)

    converged_index = next(
        i for i, line in enumerate(out_data)
        if "***        THE OPTIMIZATION HAS CONVERGED     ***" in line
    )
    conv_data = out_data[converged_index:]
    charge_index = next(i for i, line in enumerate(conv_data) if "Total Charge" in line)
    charge_str = conv_data[charge_index].split("....")[1].strip()
    mult_str   = conv_data[charge_index + 1].split("....")[1].strip()
    charge = float(charge_str)
    mult   = float(mult_str)

    # 総双極子モーメント
    di_mo_index = next(i for i, line in enumerate(conv_data) if "Total Dipole Moment" in line)
    dipole_moment = conv_data[di_mo_index].split(":")[1].strip().split()
    di_moment = ", ".join(dipole_moment)

    # ポラリザビリティ (raw cartesian tensor)
    pol_index = next(i for i, line in enumerate(out_data) if "The raw cartesian tensor (atomic units):" in line)
    pol_row = out_data[pol_index + 1 : pol_index + 4]
    pol_line = list(range(6))  # [xx, xy, xz, yy, yz, zz] のように格納
    for i, line_ in enumerate(pol_row):
        parts = line_.split()
        if i == 0:
            pol_line[0] = parts[0]  # xx
            pol_line[1] = parts[1]  # xy
            pol_line[2] = parts[2]  # xz
        elif i == 1:
            pol_line[3] = parts[1]  # yy
            pol_line[4] = parts[2]  # yz
        elif i == 2:
            pol_line[5] = parts[2]  # zz
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

    # .minfo の書き出し
    with open(minfo_path, "w", encoding="UTF-8") as f:
        f.write("[ Electronic Data ]\n")
        f.write("Energy\n")
        f.write(str(energy) + "\n")
        f.write("Charge\n")
        f.write(str(charge) + "\n")
        f.write("Multiplicity\n")
        f.write(str(mult) + "\n")
        f.write("Gradient\n")
        f.write(str(atom_num * 3) + "\n")
        f.write(grad + "\n")
        f.write("Hessian\n")
        hessian_count = atom_num*3 * (atom_num*3 + 1) // 2
        f.write(str(hessian_count) + "\n")
        f.write(upper_triangle + "\n")

        f.write("Dipole Moment\n3\n")
        f.write(di_moment + "\n")

        f.write("Polarizability\n6\n")
        f.write(pol + "\n")

        f.write("Dipole Derivative\n")
        f.write(str(atom_num*3*3) + "\n")
        f.write(dipole + "\n")

    print(f"[DONE] {work_dir}: {mol}.minfo を出力しました.")


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

    # 処理対象となるディレクトリを探索
    current_dir = os.getcwd()
    for d in sorted(os.listdir(current_dir)):
        # mkqff で始まるかどうか
        if os.path.isdir(d) and d.startswith("mkqff"):
            print(f"[INFO] 処理ディレクトリ: {d}")
            process_directory(d, atom_num, line_flag)

if __name__ == "__main__":
    main()

