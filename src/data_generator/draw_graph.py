import matplotlib
matplotlib.use('Agg')  # 不使用 GUI 的 headless 後端
import networkx as nx
import matplotlib.pyplot as plt
import re



# ====== 資料路徑 ======
input_file = "dataset/output/graph.txt"

# ====== 解析資料並建構圖 ======
G = nx.Graph()
node_type = {}  # 用來記錄每個 node 是什麼類型

with open(input_file, "r") as f:
    lines = f.readlines()

i = 1
while i < len(lines):
    line = lines[i].strip()
    
    if line.startswith("node"):
        node_id = int(line.split()[1])
        i += 1
        
        # 解析類型行
        type_line = lines[i].strip()
        types = re.findall(r'gs1 (\d+)|gs2 (\d+)|sat (\d+)', type_line)
        label = []
        for t in types:
            if t[0]:
                label.append("gs1")
            elif t[1]:
                label.append("gs2")
            elif t[2]:
                label.append("sat")
        node_type[node_id] = ",".join(label)
        G.add_node(node_id, type=node_type[node_id])
        
        i += 1  # skip weight line
        i += 1

        # 解析鄰居行
        neighbors = list(map(int, lines[i].strip().split()[1:]))
        for n in neighbors:
            G.add_edge(node_id, n)

    i += 1

# ====== 畫圖 ======
plt.figure(figsize=(16, 12))
pos = nx.spring_layout(G, seed=42)

# 不同類型用不同顏色
color_map = []
for node in G.nodes():
    t = node_type.get(node, "")
    if "gs1" in t:
        color_map.append("skyblue")
    elif "gs2" in t:
        color_map.append("orange")
    elif "sat" in t:
        color_map.append("lightgreen")
    else:
        color_map.append("gray")

nx.draw_networkx_nodes(G, pos, node_color=color_map, node_size=50)
nx.draw_networkx_edges(G, pos, alpha=0.3, width=0.5)
nx.draw_networkx_labels(G, pos, labels={n: str(n) for n in G.nodes if n < 3}, font_size=8)

plt.title("Graph Visualization", fontsize=16)
plt.axis("off")
plt.savefig("graph.png", dpi=300)
plt.close()

