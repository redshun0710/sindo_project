#!/usr/bin/env python3
import numpy as np
def add_weights_to_atomic_line(atomic_line, atom_weight_line):
    for i, (atom, weight_data) in enumerate(zip(atomic_line, atom_weight_line)):
        if atom[0] == weight_data[0]: 
            weight = weight_data[1]
            atom.insert(1, weight)
    return atomic_line

def main():
    mol = input("ディレクトリを入力してください")
    line_flag = input("直線分子かどうかTrue or False")
    atom_num = input("原子数")
    # --- hess 読み込み ---
    with open(f"{mol}.hess", "r", encoding="UTF-8") as f:
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
    # --- 最終エネルギーと荷電などを出力ファイルから取得 ---
    with open(f"{mol}.out","r",encoding="UTF-8") as f:
        data = f.readlines()
    energy_index = next(i for i, line in enumerate(data) if "*** OPTIMIZATION RUN DONE ***" in line) - 3
    energy = float(data[energy_index].split("FINAL SINGLE POINT ENERGY")[1].split("\n")[0])
    converged_index = next(i for i, line in enumerate(data) if "***        THE OPTIMIZATION HAS CONVERGED     ***" in line)
    conv_data = data[converged_index:]
    charge_index = next(i for i, line in enumerate(conv_data) if "Total Charge" in line)
    charge = float(conv_data[charge_index].split("....")[1].split("\n")[0])
    mult = float(conv_data[charge_index+1].split("....")[1].split("\n")[0])

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
if __name__ == "__main__":
    main()

