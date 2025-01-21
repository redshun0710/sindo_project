#!/usr/bin/env python3
import numpy as np
def extract_numbers(data):
    all_numbers = []
    for line in data:
        parts = line.split()
        numbers = [float(p) for p in parts if "E" in p or "." in p]
        if numbers:
            all_numbers.append(numbers)
    return all_numbers

def extract_hessian_custom_order(lines, atom_num):
    hessian_start = lines.index("$hessian\n") + 3
    N = 3 * atom_num
    
    hessian_end = hessian_start + 1 + ((atom_num**2 // 5 + 1) * atom_num**2)
    hessian_lines = lines[hessian_start:hessian_end]
    
    row_matrix = extract_numbers(hessian_lines)
    
    group_size = atom_num**2
    num_groups = len(row_matrix) // group_size
    groups = [row_matrix[i * group_size:(i + 1) * group_size] for i in range(num_groups)]
    
    matrix = []
    for i in range(group_size):
        combined_row = []
        for group in groups:
            combined_row.extend(group[i])
        matrix.append(combined_row)

    custom_order = []
    for i in range(len(matrix)):
        for j in range(i + 1):
            custom_order.append(str(matrix[i][j]))
    
    return custom_order

def add_weights_to_atomic_line(atomic_line, atom_weight_line):
    for i, (atom, weight_data) in enumerate(zip(atomic_line, atom_weight_line)):
        if atom[0] == weight_data[0]: 
            weight = weight_data[1]
            atom.insert(1, weight)
    return atomic_line

def main():
    mol = input("ディレクトリを入力してください")
    line_flag = input("直線分子かどうかTrue or False")

    # --- xyz 読み込み ---
    with open(f"{mol}_opt.xyz", "r", encoding="UTF-8") as f:
        xyzdata = f.readlines()
    atom_num = int(xyzdata[0])
    
    atomic_line = xyzdata[2 : 2 + atom_num]
    atomic_line = [line.strip().split() for line in atomic_line]

    # --- hess 読み込み ---
    with open(f"{mol}_opt.hess", "r", encoding="UTF-8") as f:
        hessdata = f.readlines()

    atom_weight_index = hessdata.index("$atoms\n")+2
    atom_weight_line = hessdata[atom_weight_index : atom_weight_index + atom_num]
    atom_weight_line = [line.strip().split() for line in atom_weight_line]
    
    atomic_line = add_weights_to_atomic_line(atomic_line, atom_weight_line)

    dipole_index = hessdata.index("$dipole_derivatives\n")+2
    dipole_line = hessdata[dipole_index : dipole_index + atom_num*3]
    dipole_line = [line.strip().split() for line in dipole_line]
    dipole_data = [
        dipole_line[row][col]
        for col in range(3)
        for row in range(atom_num*3)
    ]
    dipole = "\n".join(", ".join(dipole_data[i:i+5]) for i in range(0, len(dipole_data), 5))
    
    vib_freq_index = hessdata.index("$vibrational_frequencies\n")+2
    vib_freq_line = hessdata[vib_freq_index : vib_freq_index + atom_num*3]
    vib_freq = [line.strip().split()[1] for line in vib_freq_line]

    if line_flag == True:
        rot_num = 2
        tra = list(range(3))
        rot = list(range(rot_num))
        vib = list(range(atom_num*3-3-rot_num))
        for i, val in enumerate(vib_freq):
            if i <=2:
                tra[i] = vib_freq[i]
            elif (i > 2) and (i <= 2+rot_num):
                rot[i-3] = vib_freq[i]
            else:
                vib[i-3-rot_num] = vib_freq[i]
    else:
        rot_num = 3
        tra = list(range(3))
        rot = list(range(rot_num))
        vib = list(range(atom_num*3-3-rot_num))
        for i, val in enumerate(vib_freq):
            if i <=2:
                tra[i] = vib_freq[i]
            elif (i > 2) and (i <= 2+rot_num):
                rot[i-3] = vib_freq[i]
            else:
                vib[i-3-rot_num] = vib_freq[i]
    
    trans = ", ".join(val for val in tra)
    rotat = ", ".join(val for val in rot)
    vibra = ", ".join(val for val in vib)

    # --- normal_modes 読み込み ---
    filename = f"{mol}_opt.hess"
    mat_size = atom_num * 3
    with open(filename, 'r', encoding='UTF-8') as f:
        lines = f.readlines()
    
    start_index = None
    for i, line in enumerate(lines):
        if line.strip().startswith('$normal_modes'):
            start_index = i
            break
    
    mode_line_index = None
    for i in range(start_index + 1, len(lines)):
        parts = lines[i].split()
        if len(parts) == 2:
            try:
                r1, r2 = map(int, parts)
                if (r1 == mat_size) and (r2 == mat_size):
                    mode_line_index = i
                    break
            except:
                pass
    
    big_array = np.zeros((mat_size, mat_size))
    
    block1_header_index = mode_line_index + 1
    block1_data_start = block1_header_index + 1
    block2_header_index = block1_data_start + mat_size
    block2_data_start = block2_header_index + 1

    # --- normal_modes: ブロック1 ---
    for row in range(mat_size):
        line = lines[block1_data_start + row]
        parts = line.split()
        vals = [float(x.replace('E','e')) for x in parts[1:1+5]]
        for col in range(5):
            big_array[row, col] = vals[col]
    
    # --- normal_modes: ブロック2 ---
    for row in range(mat_size):
        line = lines[block2_data_start + row]
        parts = line.split()
        vals = [float(x.replace('E','e')) for x in parts[1:1+4]]
        for col in range(4):
            big_array[row, 5 + col] = vals[col]
    
    tra_mode = big_array[:, :3]
    rot_mode = big_array[:, 3:3 + rot_num]
    vib_mode = big_array[:, 3 + rot_num:]
    
    tra_modes = [tra_mode[:, i].tolist() for i in range(tra_mode.shape[1])]
    rot_modes = [rot_mode[:, i].tolist() for i in range(rot_mode.shape[1])]
    vib_modes = [vib_mode[:, i].tolist() for i in range(vib_mode.shape[1])]

    # --- Hessian (custom order) ---
    upper_triangle_line = extract_hessian_custom_order(hessdata, atom_num)
    upper_triangle = "\n".join(", ".join(upper_triangle_line[i:i+5]) for i in range(0, len(upper_triangle_line), 5))

    # --- 最終エネルギーと荷電などを出力ファイルから取得 ---
    with open(f"{mol}_opt.out","r",encoding="UTF-8") as f:
        data = f.readlines()
    energy_index = next(i for i, line in enumerate(data) if "*** OPTIMIZATION RUN DONE ***" in line) - 3
    energy = float(data[energy_index].split("FINAL SINGLE POINT ENERGY")[1].split("\n")[0])
    converged_index = next(i for i, line in enumerate(data) if "***        THE OPTIMIZATION HAS CONVERGED     ***" in line)
    conv_data = data[converged_index:]
    charge_index = next(i for i, line in enumerate(conv_data) if "Total Charge" in line)
    charge = float(conv_data[charge_index].split("....")[1].split("\n")[0])
    mult = float(conv_data[charge_index+1].split("....")[1].split("\n")[0])

    # --- 原子番号を追加 ---
    nuc_index = next(i for i, line in enumerate(conv_data) if "CARTESIAN COORDINATES (A.U.)" in line)
    inserted_indices = set()
    for i in range(atom_num):
        nuc_num = conv_data[nuc_index + 3 + i].split()[1:3]
        # 原子番号を int に変換したうえで再度 str() に
        nuc_num[1] = str(int(float(nuc_num[1])))
        for j, atomic in enumerate(atomic_line):
            if nuc_num[0] == atomic[0] and j not in inserted_indices:
                atomic.insert(1, nuc_num[1])
                inserted_indices.add(j)
                break
    
    atomic_xyz = [
        ", ".join(parts) for parts in atomic_line
    ]

    # --- 総双極子モーメント ---
    di_mo_index = next(i for i, line in enumerate(conv_data) if "Total Dipole Moment" in line)
    dipole_moment = conv_data[di_mo_index].split(":")[1].strip("\n").split()
    di_moment = ", ".join(val for val in dipole_moment)

    # --- ポーラリザビリティ (raw cartesian tensor) ---
    pol_index = next(i for i, line in enumerate(data) if "The raw cartesian tensor (atomic units):" in line)
    pol_row = data[pol_index+1 : pol_index+4]
    pol_line = list(range(6))
    for i,line in enumerate(pol_row):
        nline = line.split()
        if i == 0:
            pol_line[0] = str(nline[0])
            pol_line[1] = str(nline[1])
            pol_line[2] = str(nline[2])
        elif i == 1:
            pol_line[3] = str(nline[1])
            pol_line[4] = str(nline[2])
        elif i == 2:
            pol_line[5] = str(nline[2])
    pol = "\n".join(", ".join(str(x) for x in pol_line[i:i+5]) for i in range(0, len(pol_line), 5))

    # --- 勾配 (grad) ---
    with open(f"{mol}_opt.engrad","r",encoding="UTF-8") as f:
        g_data = f.readlines()
    grad_index = next(i for i, line in enumerate(g_data) if "The current gradient" in line)
    grad_line = g_data[grad_index+2:grad_index+2+atom_num*3]
    grad_line = [line.strip() for line in grad_line]
    grad = "\n".join(", ".join(grad_line[i:i+5]) for i in range(0, len(grad_line), 5))

    # --- .minfo ファイルの書き出し ---
    with open(f"{mol}.minfo", "w", encoding="UTF-8") as f:
        f.write("[ Atomic Data ]\n")
        f.write(str(atom_num) + "\n")
        for line in atomic_xyz:
            f.write(f" {line}\n")
        f.write("\n")
        f.write("[ Electronic Data ]\n")
        f.write("Energy\n")
        f.write(str(energy))
        f.write("\nCharge\n")
        f.write(str(charge))
        f.write("\nMultiplicity\n")
        f.write(str(mult))
        f.write("\nGradient\n")
        f.write(str(atom_num*3))
        f.write("\n")
        f.write(grad)
        f.write("\n")
        f.write("Hessian\n")
        f.write(str(int(atom_num**2*(atom_num**2+1)/2)))
        f.write("\n")
        f.write(upper_triangle)
        f.write("\nDipole Moment\n3\n")
        f.write(di_moment)
        f.write("\nPolarizability\n")
        f.write("6\n")
        f.write(pol)
        f.write("\n")
        f.write("Dipole Derivative\n")
        f.write(str(atom_num*3*3))
        f.write("\n")
        f.write(dipole)
        f.write("\n[ Vibrational Data ]")
        f.write("\nNormal modes")
        f.write("\nTranslational Frequency\n")
        f.write("3\n")
        f.write(trans)
        f.write("\nTranslational vector\n")
        for i, val in enumerate(tra_modes):
            f.write(f"T {i+1}\n")
            f.write(str(atom_num*3))
            f.write("\n")
            v = "\n".join(", ".join(str(x) for x in val[i:i+5]) for i in range(0, len(val), 5))
            f.write(v)
            f.write("\n")
        f.write("\nRotational Frequency\n")
        f.write(f"{rot_num}\n")
        f.write(rotat)
        f.write("\nRotational vector\n")
        for i, val in enumerate(rot_modes):
            f.write(f"R {i+1}\n")
            f.write(str(atom_num*3))
            f.write("\n")
            v = "\n".join(", ".join(str(x) for x in val[i:i+5]) for i in range(0, len(val), 5))
            f.write(v)
            f.write("\n")
        f.write("\nVibrational Frequency\n")
        f.write(f"{len(vib_modes)}")
        f.write("\n")
        f.write(vibra)
        f.write("\nVibrational vector\n")
        for i, val in enumerate(vib_modes):
            f.write(f"Mode {i+1}\n")
            f.write(str(atom_num*3))
            f.write("\n")
            v = "\n".join(", ".join(str(x) for x in val[i:i+5]) for i in range(0, len(val), 5))
            f.write(v)
            f.write("\n")

if __name__ == "__main__":
    main()

