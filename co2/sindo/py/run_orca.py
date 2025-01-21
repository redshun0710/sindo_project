#!/usr/bin/env python3

import os
import subprocess

def main():
    xyzfile = input("xyzファイルのディレクトリ")
    # 例: orca 実行コマンド (パスを指定する場合は絶対パスで書く)
    orca_cmd = "orca"

    # ========================================
    # 1. xyz ファイルを読み込み、ブロックに分解
    # ========================================
    with open(xyzfile, "r", encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f]

    idx = 0
    total_lines = len(lines)

    while idx < total_lines:
        # --- (a) 原子数を取得 ---
        natom_line = lines[idx].strip()
        # 空行などをスキップする場合は適宜条件分岐を追加
        natom = int(natom_line)
        idx += 1

        # --- (b) ブロック名を取得 (mkqff-eq など) ---
        block_name = lines[idx].strip()
        idx += 1

        # --- (c) 原子座標を取得 ---
        coord_lines = lines[idx : idx + natom]
        idx += natom

        # ========================================
        # 2. ディレクトリ作成 & ORCA入力ファイル作成
        # ========================================
        # ディレクトリ名はブロック名に対応
        work_dir = block_name
        os.makedirs(work_dir, exist_ok=True)

        # --- ORCA 入力ファイルを準備 ---
        # ここでは例として B3LYP/def2-SVP で Frequency 計算を行い、
        # 与えられた座標を固定したままエネルギー & Hessian を求める想定
        orca_input = make_orca_input(block_name, coord_lines)

        inp_path = os.path.join(work_dir, "job.inp")
        with open(inp_path, "w", encoding="utf-8") as f_inp:
            f_inp.write(orca_input)

        # ========================================
        # 3. ORCA で計算実行
        # ========================================
        # 計算コマンド例:  orca job.inp > job.out
        # subprocess.run の引数でカレントディレクトリを指定
        cmd = f"{orca_cmd} job.inp > job.out"
        print(f"[INFO] Running ORCA in {work_dir} ...")
        subprocess.run(cmd, shell=True, cwd=work_dir)

        # 以降、必要に応じて .out や .hess, .engrad などを読み取る処理を書いてもよい

    print("All jobs submitted/finished.")


def make_orca_input(block_name, coord_lines):
    charge = 0
    multiplicity = 1
    geom_block = "\n".join(coord_lines)
    orca_input = f"""! B3LYP def2-SVP TightSCF Freq
%output
  PrintLevel 5
end
* xyz {charge} {multiplicity}
{geom_block}
*
%elprop
 Polar 1
end
"""
    return orca_input


if __name__ == "__main__":
    main()

