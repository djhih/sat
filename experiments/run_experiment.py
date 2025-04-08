import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Microsoft JhengHei", "SimHei", "Arial Unicode MS", "sans-serif"]
matplotlib.rcParams["axes.unicode_minus"] = False

def run_command(cmd, description="", timeout=None):
    print(description)
    try:
        subprocess.run(cmd, check=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        print(f"執行命令 {cmd} 超過 {timeout} 秒，強制停止。")
    except subprocess.CalledProcessError as e:
        print(f"執行命令 {cmd} 時發生錯誤：{e}")
        raise

def run_data_generator(num_sat, num_req, sat_filename, gs_filename, output_filename):
    """
    編譯並執行 C++ 版的資料產生器：
      - 先利用 Python 產生衛星與地面站資料
      - 再編譯並執行 C++ 資料產生器
    """
    # 生成 n 個衛星資料
    run_command(
        ["python3", "dataset/code/tle_gen_n.py", str(num_sat)],
        f"產生 {num_sat} 筆衛星資料..."
    )
    
    # 生成地面站資料
    run_command(
        ["python3", "dataset/code/USA_gov_location.py"],
        "產生地面站資料..."
    )
    
    # 定義資料產生器的原始碼與執行檔路徑
    data_gen_src = Path("src") / "data_generator" / "gen_v2.cpp"
    data_gen_exe = Path("bin") / "gen_v2.exe"
    
    run_command(
        ["g++", str(data_gen_src), "-o", str(data_gen_exe)],
        f"編譯 {data_gen_src}..."
    )
    
    run_command(
        [str(data_gen_exe), str(num_req), sat_filename, gs_filename, output_filename],
        "執行資料產生器..."
    )
    print(f"資料已寫入 {output_filename}\n")

def compile_and_run_algorithms(dataset_file, dataset_id):
    """
    編譯並執行 C++ 演算法：
      - 編譯 src/algorithms/ 下的 C++ 檔案
      - 執行後將結果依據 dataset_id 寫入 dataset/output/ 下
      - 若某演算法執行超過 2 分鐘 120 秒，則強制停止
    """
    algorithms = {
        "greedy": (Path("src") / "algorithms" / "greedy.cpp",
                   Path("dataset") / "output" / f"res_greedy_{dataset_id}.txt"),
        "heuristic": (Path("src") / "algorithms" / "heuristic.cpp",
                      Path("dataset") / "output" / f"res_he_{dataset_id}.txt"),
        "max_independent": (Path("src") / "algorithms" / "max_independent_set.cpp",
                            Path("dataset") / "output" / f"res_max_{dataset_id}.txt"),
        "imp": (Path("src") / "algorithms" / "imp.cpp",
                Path("dataset") / "output" / f"res_imp_{dataset_id}.txt")
    }

    for algo, (src_path, res_file) in algorithms.items():
        exe_path = Path("bin") / f"{algo}.exe"
        run_command(
            ["g++", str(src_path), "-o", str(exe_path), "-std=c++17"],
            f"編譯 {src_path}..."
        )
        
        # 若演算法執行超過 2 分鐘則會超時
        run_command(
            [str(exe_path), dataset_file, str(res_file)],
            f"執行 {exe_path}...",
        )

    return algorithms

def read_all_results(dataset_ids):
    """
    讀取各 dataset 的演算法結果檔（只處理以 "accept" 開頭的行），
    解析後回傳一個巢狀字典，其結構如下：
      { dataset_id: { "greedy": [...], "heuristic": [...], "max_independent": [...], "imp": [...] } }
    """
    results = {}
    algo_suffixes = {
        "greedy": "greedy",
        "heuristic": "he",
        "max_independent": "max",
        "imp": "imp"
    }
    
    for ds_id in dataset_ids:
        results[ds_id] = {}
        for algo, suffix in algo_suffixes.items():
            res_file = Path("dataset") / "output" / f"res_{suffix}_{ds_id}.txt"
            print(f"讀取結果檔：{res_file}")
            data = []
            try:
                with open(res_file, "r") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("accept"):
                            tokens = line.split()
                            try:
                                parsed = [
                                    float(tokens[2]),
                                    float(tokens[4]),
                                    float(tokens[6]),
                                    float(tokens[9])
                                ]
                                data.append(parsed)
                            except (IndexError, ValueError) as e:
                                print(f"解析行失敗: '{line}', error: {e}")
                                continue
            except FileNotFoundError:
                print(f"找不到檔案: {res_file}")
            results[ds_id][algo] = data
    return results

def plot_grouped_total_rates_by_algorithm(dataset_ids):
    """
    讀取各 dataset 中 "Total generation rate:" 數值，
    並以分組柱狀圖呈現不同演算法的比較：
      - x 軸為 dataset_id，每個 dataset 依序有四個柱子
    """
    algo_suffixes = {
        "greedy": "greedy",
        "heuristic": "he",
        "max_independent": "max",
        "imp": "imp"
    }
    
    total_rates = { algo: [] for algo in algo_suffixes }
    
    for ds_id in dataset_ids:
        for algo, suffix in algo_suffixes.items():
            res_file = Path("dataset") / "output" / f"res_{suffix}_{ds_id}.txt"
            tot_rate = None
            try:
                with open(res_file, "r") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("Total generation rate:"):
                            parts = line.split(":")
                            if len(parts) >= 2:
                                tot_rate = float(parts[1].strip())
                            break
            except Exception as e:
                print(f"讀取檔案 {res_file} 失敗: {e}")
            if tot_rate is None:
                print(f"檔案 {res_file} 中未找到 'Total generation rate:'，預設 0.0")
                tot_rate = 0.0
            total_rates[algo].append(tot_rate)
    
    # 將 dataset_ids 轉換為整數（假設 dataset_ids 為 "1", "2", ...）
    x = np.array([int(ds) for ds in dataset_ids])
    
    n_algs = len(algo_suffixes)
    width = 0.8 / n_algs
    offsets = np.arange(n_algs) * width - (0.8 - width) / 2
    
    fig, ax = plt.subplots()
    algo_names = list(algo_suffixes.keys())
    for i, algo in enumerate(algo_names):
        positions = x + offsets[i]
        ax.bar(positions, total_rates[algo], width=width, label=algo)
    
    ax.set_xlabel("Dataset ID")
    ax.set_ylabel("Total Generation Rate")
    ax.set_title("各演算法之 Total Generation Rate 比較 (分組柱狀圖)")
    ax.set_xticks(x)
    ax.set_xticklabels(dataset_ids)
    ax.legend()
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)
    
    plt.tight_layout()
    filename = "grouped_total_generation_rate_by_algorithm.png"
    plt.savefig(filename)
    print(f"圖表已儲存：{filename}")
    plt.show()

def main():
    # 檔案路徑設定
    sat_all_filename = "dataset/code/output/satellite_coordinates_all.txt"
    sat_n_filename = "dataset/code/output/satellite_coordinates_selected.txt"
    op_filename = "dataset/raw/dataset"  # 後面會接上 .txt
    gs_filename = "dataset/code/output/gs_loc.txt"

    dataset_ids = []
    
    # 產生所有衛星資料
    # run_command(
    #     ["python3", "dataset/code/tle_gen_all.py"],
    #     "產生所有衛星資料..."
    # )
    
    # 執行多組實驗
    for idx in range(1, 5):
        output_filename = op_filename + str(idx) + ".txt"
        num_req = 100
        num_sat = idx * 5
        
        dataset_ids.append(str(idx))
        print(f"\n===== 執行實驗：{idx} =====")
        run_data_generator(num_sat, num_req, sat_n_filename, gs_filename, output_filename)
        compile_and_run_algorithms(output_filename, idx)
    
    # 讀取各 dataset 的結果（若需要進一步處理可使用）
    results = read_all_results(dataset_ids)
    # 根據 dataset_ids 畫出分組柱狀圖
    plot_grouped_total_rates_by_algorithm(dataset_ids)

if __name__ == "__main__":
    main()
