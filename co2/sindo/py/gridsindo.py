import os
import subprocess

def parse_xyz(filename):
    """
    makeGrid.xyz のように
        3
        label
        C x y z
        O x y z
        O x y z
        3
        label
        ...
    のブロックが続くファイルをパースして、
    [(label, [(atom, x, y, z), ... ]), ... ] のリストを返す。
    """
    data = []
    with open(filename, 'r') as f:
        lines = f.read().strip().splitlines()
    
    i = 0
    while i < len(lines):
        # 原子数
        natoms = int(lines[i].strip())
        i += 1
        
        # ラベル(例: mkg-eq, mkg-q1-11-0 など)
        label = lines[i].strip()
        i += 1
        
        # 座標を読み込む
        coords = []
        for _ in range(natoms):
            parts = lines[i].split()
            atom = parts[0]
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
            coords.append((atom, x, y, z))
            i += 1
        
        data.append((label, coords))
    return data


def write_orca_input(filepath, coords):
    """
    ORCAの入力ファイルを簡単に書き出す例。
    必要に応じてキーワードや理論レベルを調整してください。
    """
    with open(filepath, 'w') as f:
        # 計算レベル等は適宜変更してください
        f.write("! B3LYP def2-SVP TightSCF\n\n")
        f.write("* xyz 0 1\n")
        for atom, x, y, z in coords:
            f.write(f" {atom:2s}  {x:15.10f}  {y:15.10f}  {z:15.10f}\n")
        f.write("*\n")


def main():
    xyzfile = "../makeGrid.xyz"
    xyz_data = parse_xyz(xyzfile)
    
    for label, coords in xyz_data:
        if label == "mkg-eq":
            # eq の場合は特別に eq/step0/ とする
            vib_name = "eq"
            step_name = "step0"
        else:
            # 例: label = "mkg-q1-11-0"
            # "-" でスプリットすると ["mkg", "q1", "11", "0"] のようになる想定
            parts = label.split("-")
            # parts[1]: "q1", "q2", ...
            vib_name = parts[1]  # q1, q2, q3, q4 など
            # parts[3]: ステップ番号 "0", "1", ...
            step_name = "step" + parts[3]
        
        # ディレクトリ作成 (例: "q1/step0" など)
        outdir = os.path.join(vib_name, step_name)
        os.makedirs(outdir, exist_ok=True)
        
        # ORCA 入力ファイル
        inp_file = os.path.join(outdir, "calc.inp")
        # ORCA 出力ファイル
        out_file = os.path.join(outdir, "calc.out")
        
        # 入力ファイルを書き出し
        write_orca_input(inp_file, coords)
        
        # ORCA を実行 (計算させる場合)
        with open(out_file, 'w') as outf:
            # ここで ORCA を呼び出す。パスが必要なら "path/to/orca" に書き換え。
            subprocess.run(["orca", inp_file], stdout=outf, stderr=outf)


if __name__ == "__main__":
    main()

