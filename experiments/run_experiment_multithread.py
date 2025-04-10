import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_data_generator(num_rings, num_sat, gs_filename, output_filename):
    """
    編譯並執行 C++ 版的資料產生器：
      - 原始碼位於 src/data_generator/cpp_data_generator.cpp
      - 編譯後的執行檔存放在 bin/data_generator.exe
      - argv[1] 是 number of rings
      - argv[2] 是 number of sat each rings
      - argv[3] 是 gs_filename (產生的 dataset 檔名)
      - argv[4] 是 output filename
    """
    data_gen_src = os.path.join("src", "data_generator", "cpp_data_generator.cpp")
    data_gen_exe = os.path.join("bin", "data_generator.exe")
    
    print(f"編譯資料產生器：{data_gen_src} -> {data_gen_exe}")
    subprocess.run(["g++", data_gen_src, "-o", data_gen_exe], check=True)
    
    print(f"執行資料產生器：{data_gen_exe} {num_rings} {num_sat} {gs_filename} {output_filename}")
    subprocess.run([data_gen_exe, str(num_rings), str(num_sat), gs_filename, output_filename], check=True)

def compile_and_run_algorithms(dataset_file, dataset_id):
    """
    編譯並執行 C++ 演算法：
      - 原始碼位於 src/algorithms/
      - 編譯後的執行檔放在 bin/ 中
      - 執行時傳入 dataset filename (argv[1])
      - 演算法的標準輸出將依據 dataset_id 寫入 dataset/output/ 下不同檔案
    """
    algorithms = {
        "greedy": (os.path.join("src", "algorithms", "greedy.cpp"),
                   os.path.join("dataset", "output", f"res_greedy_{dataset_id}.txt")),
        "heuristic": (os.path.join("src", "algorithms", "heuristic.cpp"),
                      os.path.join("dataset", "output", f"res_he_{dataset_id}.txt")),
        "max_independent": (os.path.join("src", "algorithms", "max_independent_set.cpp"),
                            os.path.join("dataset", "output", f"res_max_{dataset_id}.txt")),
        "imp": (os.path.join("src", "algorithms", "imp_multithread.cpp"),
                os.path.join("dataset", "output", f"res_imp_{dataset_id}.txt")),
    }

    for algo, (src_path, res_file) in algorithms.items():
        exe_path = os.path.join("bin", algo + ".exe")
        print(f"編譯 {src_path} 成 {exe_path} ...")
        subprocess.run(["g++", src_path, "-o", exe_path], check=True)
        
        print(f"執行 {exe_path}，讀取 dataset: {dataset_file}，結果將寫入 {res_file} ...")
        subprocess.run([exe_path, dataset_file, res_file], text=True, check=True)

    return algorithms

def read_all_results(dataset_ids):
    """
    讀取不同 dataset 的所有演算法結果檔，檔案格式範例如下：
      accept gs1 1 gs2 2 sat 2 gen rate 40
      accept gs1 0 gs2 3 sat 1 gen rate 24
      Total generation rate: 64

    本函數只處理以 "accept" 開頭的行，解析後回傳 [gs1, gs2, sat, gen_rate] 四個數值，
    並將結果存入一個巢狀字典，結構如下：
      { dataset_id: { "greedy": [...], "heuristic": [...], "max_independent": [...], "imp": [...] } }
    """
    results = {}
    algo_suffixes = {
        "greedy": "greedy",
        "heuristic": "he",
        "max_independent": "max",
        "imp": "imp"
    }
    
    for dataset_id in dataset_ids:
        results[dataset_id] = {}
        for algo, suffix in algo_suffixes.items():
            res_file = os.path.join("dataset", "output", f"res_{suffix}_{dataset_id}.txt")
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
            results[dataset_id][algo] = data
    return results

def plot_grouped_total_rates_by_algorithm(dataset_ids):
    """
    針對每個 dataset（以 dataset_id 為識別字），
    讀取四個演算法結果檔中 "Total generation rate:" 行的數值，
    並以分組柱狀圖呈現：
      - x 軸為 dataset_id
      - 每個 dataset 有四個柱子，分別代表不同演算法
    """
    algo_suffixes = {
        "greedy": "greedy",
        "heuristic": "he",
        "max_independent": "max",
        "imp": "imp"
    }
    
    total_rates = { algo: [] for algo in algo_suffixes }
    
    for ds in dataset_ids:
        for algo, suffix in algo_suffixes.items():
            res_file = os.path.join("dataset", "output", f"res_{suffix}_{ds}.txt")
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
                print(f"檔案 {res_file} 中未找到 Total generation rate，預設 0.0")
                tot_rate = 0.0
            total_rates[algo].append(tot_rate)
    
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
    ax.set_title("Total Generation Rate Comparison by Algorithm (Grouped Bar Chart)")
    ax.set_xticks(x)
    ax.set_xticklabels(dataset_ids)
    ax.legend()
    ax.grid(True, axis="y", linestyle="--", alpha=0.7)
    
    plt.tight_layout()
    filename = "grouped_total_generation_rate_by_algorithm.png"
    plt.savefig(filename)
    print(f"圖表已儲存：{filename}")
    plt.show()

def process_dataset(params, dataset_id):
    """
    單一 dataset 的完整處理流程：
      1. 產生 dataset（呼叫 run_data_generator）
      2. 編譯並執行各演算法（呼叫 compile_and_run_algorithms）
    傳入 dataset_id 用作輸出檔名的識別字。
    """
    num_rings = params["num_rings"]
    num_sat = params["num_sat"]
    gs_filename = params["gs_filename"]
    output_filename = params["op_filename"]
    print(f"\n===== 執行實驗：{dataset_id} =====")
    run_data_generator(num_rings, num_sat, gs_filename, output_filename)
    compile_and_run_algorithms(output_filename, dataset_id)
    return dataset_id

if __name__ == "__main__":
    # 定義多組 dataset 的參數組合
    datasets_params = [
        {"num_rings": 15, "num_sat": 10, "gs_filename": "dataset/code/output/gs_loc.txt", "op_filename": "dataset/raw/dataset1.txt"},
        {"num_rings": 15, "num_sat": 15, "gs_filename": "dataset/code/output/gs_loc.txt", "op_filename": "dataset/raw/dataset2.txt"},
        {"num_rings": 20, "num_sat": 20, "gs_filename": "dataset/code/output/gs_loc.txt", "op_filename": "dataset/raw/dataset3.txt"},
        {"num_rings": 25, "num_sat": 25, "gs_filename": "dataset/code/output/gs_loc.txt", "op_filename": "dataset/raw/dataset4.txt"},
        {"num_rings": 30, "num_sat": 25, "gs_filename": "dataset/code/output/gs_loc.txt", "op_filename": "dataset/raw/dataset5.txt"},
        {"num_rings": 30, "num_sat": 30, "gs_filename": "dataset/code/output/gs_loc.txt", "op_filename": "dataset/raw/dataset6.txt"},
    ]
    
    dataset_ids = []
    # 使用 ThreadPoolExecutor 同時處理多個 dataset
    with ThreadPoolExecutor(max_workers=len(datasets_params)) as executor:
        futures = []
        for idx, params in enumerate(datasets_params, start=1):
            futures.append(executor.submit(process_dataset, params, str(idx)))
        for future in as_completed(futures):
            ds_id = future.result()
            dataset_ids.append(ds_id)
    
    # 排序 dataset_ids（若順序有變化）
    dataset_ids.sort(key=lambda x: int(x))
    
    results = read_all_results(dataset_ids)
    plot_grouped_total_rates_by_algorithm(dataset_ids)
