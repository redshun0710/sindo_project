#!/usr/bin/env python3
import os
import re
import math
from glob import glob

def main():
    """
    minfoファイルから基準振動モードのベクトルを取得し、
    eq/calc.out（平衡構造の計算結果）と
    q1, q2, ... 各フォルダの step*/calc.out を比較して
    質量重み付き正規変位量 Q, 差分エネルギー, 差分双極子モーメントを計算。

    最終的に ../grid/q{n}.pot, ../grid/q{n}.dipole にまとめる。
    """
    # 1) ユーザー入力
    natoms_str = input("原子数を入力してください (整数): ")
    natoms = int(natoms_str.strip())

    linear_str = input("直線分子ですか？(True/False): ").strip().lower()
    is_linear = (linear_str == "true")

    mol_name = input("分子名(例: my_molecule)を入力してください(拡張子除く): ").strip()
    minfo_path = f"../{mol_name}.minfo"

    if not os.path.exists(minfo_path):
        print(f"[ERROR] {minfo_path} が見つかりません。中断します。")
        return

    # 2) minfoファイルを解析し、原子情報 & 振動モードを取得
    atomic_data, vibr_freqs, vibr_modes = parse_minfo(minfo_path)
    # atomic_data: [{'element':..., 'mass':..., 'eqx':..., 'eqy':..., 'eqz':...}, ...]
    # vibr_freqs: [652.0658..., 652.0659..., 1387.77..., 2469.87...]
    # vibr_modes: { mode_idx(1-based): [(vx1,vy1,vz1), (vx2,vy2,vz2), ...], ... }

    # 分子全体の振動モード数を計算
    if is_linear:
        nmodes = 3 * natoms - 5
    else:
        nmodes = 3 * natoms - 6

    print(f"[INFO] 予想される振動モード数: {nmodes}")
    # ただし minfo上で "Vibrational Frequency 4" のように、実際の出力数と合わない例もあるため注意

    # 3) eq/calc.out からエネルギー・双極子・計算手法(!行)・平衡座標を取得
    eq_out_path = os.path.join("eq/step0", "calc.out")
    if not os.path.exists(eq_out_path):
        print(f"[ERROR] {eq_out_path} がありません。処理を中断します。")
        return

    eq_energy, eq_dipole = parse_energy_and_dipole(eq_out_path)
    if eq_energy is None or eq_dipole is None:
        print("[ERROR] eq/step0/calc.out からエネルギー/双極子モーメント取得失敗。")
        return

    eq_basis_method, eq_coords = parse_input_section(eq_out_path)
    if len(eq_coords) != natoms:
        print(f"[WARN] eq座標の原子数({len(eq_coords)})が入力の値({natoms})と合いません。")

    # ../grid ディレクトリを作成
    grid_dir = "../grid"
    os.makedirs(grid_dir, exist_ok=True)

    # 4) 振動モードごと(q1,q2,...) に各 step*/calc.out を走査
    #    → 差分エネルギー, 差分双極子, 質量重み付き正規変位量 Q を計算→ q{n}.pot, q{n}.dipole
    for mode_index in range(1, nmodes + 1):
        q_label = f"q{mode_index}"

        # minfo上の "Mode 1", "Mode 2",... に対応
        mode_vec = vibr_modes.get(mode_index, None)
        if mode_vec is None:
            print(f"[WARN] minfo から {q_label} (Mode {mode_index}) のベクトルが見つからないためスキップ。")
            continue

        q_dir = os.path.join(q_label)
        if not os.path.isdir(q_dir):
            print(f"[WARN] ディレクトリ {q_dir} が見つかりません。スキップします。")
            continue

        step_dirs = sorted(glob(os.path.join(q_dir, "step*")))

        step_data = []  # [(step_i, Q, E_diff, (dx,dy,dz), basis_method), ...]

        for sd in step_dirs:
            step_name = os.path.basename(sd)  # "step0", "step1", ...
            m = re.search(r"step(-?\d+)", step_name)
            if m:
                step_i = int(m.group(1))
            else:
                step_i = None

            out_file = os.path.join(sd, "calc.out")
            if not os.path.exists(out_file):
                continue

            step_energy, step_dipole = parse_energy_and_dipole(out_file)
            if step_energy is None or step_dipole is None:
                continue
            e_diff = step_energy - eq_energy
            d_diff = [step_dipole[k] - eq_dipole[k] for k in range(3)]

            basis_method, step_coords = parse_input_section(out_file)
            if len(step_coords) != natoms:
                print(f"[WARN] step座標 原子数が合いません: {sd}")
                continue

            # 質量重み付き振動座標Qを計算
            Q_val = compute_weighted_displacement(
                eq_coords, step_coords, mode_vec, atomic_data
            )

            step_data.append((step_i, Q_val, e_diff, d_diff, basis_method))

        # stepインデックスで昇順ソート
        step_data.sort(key=lambda x: x[0] if x[0] is not None else 999999)

        if len(step_data) == 0:
            print(f"[WARN] {q_label} に step* データがありません。")
            continue

        final_basis_method = step_data[0][4]

        # --- q{n}.pot 出力 ---
        pot_file = os.path.join(grid_dir, f"{q_label}.pot")
        with open(pot_file, "w", encoding="utf-8") as fpot:
            fpot.write(f"{final_basis_method}\n")
            fpot.write("# Number of grids and data\n")
            fpot.write(f"    {len(step_data)}      1\n")
            fpot.write(f"#   {q_label}              Energy\n")
            for (st_i, Q_, E_, _, _) in step_data:
                fpot.write(f"    {Q_:12.8f}     {E_:15.8e}\n")

        # --- q{n}.dipole 出力 ---
        dip_file = os.path.join(grid_dir, f"{q_label}.dipole")
        with open(dip_file, "w", encoding="utf-8") as fdip:
            fdip.write(f"{final_basis_method}\n")
            fdip.write("# Number of grids and data\n")
            fdip.write(f"    {len(step_data)}      3\n")
            fdip.write(f"#   {q_label}              X                   Y                   Z\n")
            for (st_i, Q_, _, d_diff, _) in step_data:
                dx, dy, dz = d_diff
                fdip.write(f"    {Q_:12.8f}     {dx:15.8e}    {dy:15.8e}    {dz:15.8e}\n")

        print(f"[INFO] {q_label}.pot, {q_label}.dipole を出力しました。")

    print("[INFO] 全処理が完了しました。")


# -------------------------------------------------------------------
# 以下、必要な関数群
# -------------------------------------------------------------------

def parse_minfo(minfo_file):
    """
    与えられた minfoファイルから:
      - Atomic Data (mass, eqX, eqY, eqZ)
      - [ Vibrational Data ] セクション中の Normal modes, Vibrational Frequency, Vibrational vector
        を取得して返す。

    戻り値:
      atomic_data: list of dict = [
        {
          'element': str,
          'mass': float (amu),
          'eqx': float,
          'eqy': float,
          'eqz': float,
        }, ...
      ]
      vibr_freqs: list of float  (実振動数 [cm^-1] のみまたは全モード)
      vibr_modes: dict {
         mode_index(1-based): [ (vx1,vy1,vz1), (vx2,vy2,vz2), ... ]
      }
    """

    with open(minfo_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    atomic_data = []
    vibr_freqs = []
    vibr_modes = {}

    # セクションを探すフラグ
    in_atomic_data = False
    in_vibr_data = False

    # モード数などの一時保持
    nmodes_found = 0  # Vibrational Frequency で報告される数
    freq_list = []    # vibrational frequencies
    current_mode_index = None

    i = 0
    nlines = len(lines)
    while i < nlines:
        line = lines[i].rstrip("\n")

        # [ Atomic Data ] セクション検出
        if line.startswith("[ Atomic Data ]"):
            in_atomic_data = True
            in_vibr_data = False
            i += 1
            continue

        # [ Vibrational Data ] セクション検出
        if line.startswith("[ Vibrational Data ]"):
            in_atomic_data = False
            in_vibr_data = True
            i += 1
            continue

        # Atomic Data を読む
        if in_atomic_data:
            if re.match(r"^\d+$", line.strip()):
                # 原子数行(例: "3") -> その次の行から実データ
                # ただしこの数字自体は使わなくてもOK
                pass
            else:
                # 例: " C, 6, 12.01100, -0.00000264818086, 0.00000003629956, 0.00000000288260"
                parts = line.split(",")
                if len(parts) >= 6:
                    elem = parts[0].strip()
                    # parts[1]: 原子番号(6など)
                    mass = float(parts[2])
                    eqx = float(parts[3])
                    eqy = float(parts[4])
                    eqz = float(parts[5])
                    atomic_data.append({
                        "element": elem,
                        "mass": mass,
                        "eqx": eqx,
                        "eqy": eqy,
                        "eqz": eqz
                    })
            i += 1
            continue

        # Vibrational Data を読む
        if in_vibr_data:
            # "Vibrational Frequency" の行を探す
            if "Vibrational Frequency" in line:
                # 次の行に モード数(例: "4") がある -> さらに次の行に振動数リスト
                i += 1
                freq_count_line = lines[i].strip()
                nmodes_found = int(freq_count_line)
                i += 1

                # 次の行に（カンマ区切りなどで）振動数が出てくる
                freq_values_line = lines[i].strip()
                # 例: "652.0658138343734436, 652.0659170981240322, 1387.7725340583306206, 2469.8717939..."
                freq_vals = re.split(r"[,\s]+", freq_values_line)
                freq_list = [float(v) for v in freq_vals if v]  # 0-length除外
                i += 1

            # "Vibrational vector" に来るまで飛ばす
            if "Vibrational vector" in line:
                # ここから "Mode 1", "Mode 2", ... のブロックを読み取る
                i += 1
                while i < nlines:
                    vec_line = lines[i].strip()
                    if vec_line.startswith("Mode"):
                        # 例: "Mode 1"
                        parts = vec_line.split()
                        current_mode_index = int(parts[1])
                        i += 1

                        # 次の行に "9" (例: 原子数*NDOF=9)
                        nm_line = lines[i].strip()
                        nm_int = int(nm_line)  # 9
                        i += 1

                        # 次の行に 9 個の浮動小数 (x1, y1, z1, x2, y2, z2, x3, y3, z3)
                        vector_line = lines[i].strip()
                        i += 1
                        arr = re.split(r"[,\s]+", vector_line)
                        arr_floats = [float(x) for x in arr if x]

                        # (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) ... の形にまとめる
                        mode_xyz = []
                        for idx in range(0, nm_int, 3):
                            mode_xyz.append((arr_floats[idx], arr_floats[idx+1], arr_floats[idx+2]))

                        vibr_modes[current_mode_index] = mode_xyz

                    else:
                        # "Mode" 行でなければ脱出か次セクション
                        if vec_line.startswith("[") or vec_line.startswith("$") or vec_line == "":
                            # 終了
                            break
                        i += 1
                # ループ終わったら continue
            i += 1
            continue

        # どのセクションでもない行
        i += 1

    # 結果まとめ
    vibr_freqs = freq_list  # 例: [652.0, 652.0, 1387.7, 2469.8]
    return atomic_data, vibr_freqs, vibr_modes


def parse_energy_and_dipole(out_file):
    """
    ORCAのcalc.outからTotal SCF Energy (Eh) と Total Dipole Moment を抽出。
    戻り値: (energy_in_hartree, [dx,dy,dz]) or (None,None)
    """
    import re

    re_energy = re.compile(r"Total Energy\s*:\s*([-\d\.Ee+]+)\s*Eh")
    re_dip    = re.compile(r"Total Dipole Moment\s*:\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)")

    energy = None
    dipole = None
    with open(out_file, "r", encoding="utf-8") as f:
        for line in f:
            m_e = re_energy.search(line)
            if m_e:
                energy = float(m_e.group(1))
            m_d = re_dip.search(line)
            if m_d:
                dipole = [float(m_d.group(1)), float(m_d.group(2)), float(m_d.group(3))]

    if energy is None or dipole is None:
        return (None, None)
    return (energy, dipole)


def parse_input_section(out_file):
    """
    ORCAの outファイルから、INPUT FILEセクションを探し、
    - "! CCSD(T)/cc-pVTZ" のような計算手法/基底関数行(basis_method)
    - "* xyz q m" 以下にある各原子の座標
    を抽出して返す

    戻り値: (basis_method_str, [(atom, x, y, z), ...])
    """
    import re

    pattern_input_start = re.compile(r"={5,}.*INPUT FILE.*={5,}")
    pattern_input_end   = re.compile(r"\*\*\*\*END OF INPUT\*\*\*\*")

    basis_method = ""
    coords = []

    in_input_block = False

    with open(out_file, "r", encoding="utf-8") as f:
        for line in f:
            if pattern_input_start.search(line):
                in_input_block = True
                continue
            if pattern_input_end.search(line):
                in_input_block = False
                break

            if in_input_block and line.strip().startswith("|"):
                # "|  1> ! B3LYP def2-SVP" のような行をパース
                text_part = line.split(">", 1)
                if len(text_part) < 2:
                    continue
                text = text_part[1].strip()

                if text.startswith("!"):
                    # "! CCSD(T)/cc-pVTZ" など
                    basis_method = text

                elif text and not text.startswith("*"):
                    # "C  -0.0001  0.0102  0.2123" など
                    arr = text.split()
                    if len(arr) == 4:
                        atom_label = arr[0]
                        try:
                            x = float(arr[1])
                            y = float(arr[2])
                            z = float(arr[3])
                            coords.append((atom_label, x, y, z))
                        except:
                            pass

    return basis_method, coords


def compute_weighted_displacement(eq_coords, step_coords, mode_vec, atomic_data):
    """
    質量重み付き正規変位量Qを計算 (非正規化モードの場合の定義)
    Q = [ Σ_i( sqrt(m_i)* (r_i - r_i^0)·(e_i) ) ] / sqrt[ Σ_j( m_j * |e_j|^2 ) ]

    eq_coords   = [(atom_label, x0, y0, z0), ...] (N個)
    step_coords = [(atom_label, x1, y1, z1), ...] (N個)
    mode_vec    = [(vx, vy, vz), ...] (N個)  <-- minfoに書かれているMode iのベクトル
    atomic_data = [{'mass':..., 'eqx':..., 'eqy':..., 'eqz':...}, ...] (N個)
    """
    if len(eq_coords) != len(step_coords) or len(mode_vec) != len(atomic_data):
        return 0.0  # サイズ不一致時は0など

    # 分母
    denom = 0.0
    for i in range(len(mode_vec)):
        m_i = atomic_data[i]["mass"]
        vx, vy, vz = mode_vec[i]
        denom += m_i * (vx**2 + vy**2 + vz**2)
    denom = math.sqrt(denom) if denom > 1e-15 else 1e-15

    # 分子
    numerator = 0.0
    for i in range(len(eq_coords)):
        m_i = atomic_data[i]["mass"]
        vx, vy, vz = mode_vec[i]

        _, x0, y0, z0 = eq_coords[i]
        _, x1, y1, z1 = step_coords[i]
        dx = x1 - x0
        dy = y1 - y0
        dz = z1 - z0

        numerator += math.sqrt(m_i) * (dx*vx + dy*vy + dz*vz)

    return (numerator / denom)


if __name__ == "__main__":
    main()

