#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <fstream>
#include <cassert>
#include <unordered_map>
#include <algorithm>
using namespace std;

typedef pair<int, int> pii;

// 自訂 hash 結構，用來作為 unordered_map 的鍵
struct PairHash {
    size_t operator()(const pii &p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};

// 全域變數，S: 衛星數量, G: 地面站數量, R: 需求數量
int S, G, R;
// 僅保留衛星與地面站間的必要參數：訊號品質 (fidelity) 與生成率 (generation rate)
unordered_map<pii, double, PairHash> fid_gs_sat;   // (gs, sat) -> fidelity
unordered_map<pii, double, PairHash> rate_gs_sat;  // (gs, sat) -> generation rate

double fidelity_threshold = 0.5;  // 門檻值

// 衛星結構，只保留用來檢查每個 ground station 的訊號品質
struct Satellite {
    int id;
    set<int> gs_serve;   // 可服務的地面站集合（根據 fidelity 判斷）
    set<int> gsp_serve;  // 可服務的需求（地面站對）的編號集合

    // 檢查某地面站是否符合 fidelity 門檻
    bool check_fid(int gs) {
        return fid_gs_sat[{gs, id}] >= fidelity_threshold;
    }
    // 同時檢查一對地面站是否都符合門檻
    bool check_both_gs(int gs1, int gs2) {
        return check_fid(gs1) && check_fid(gs2);
    }
};

// 需求結構，只需要記錄兩個地面站的編號，以及各衛星對該需求的生成率
struct Requirement {
    int id;
    int gs1, gs2;
    int served_by_sat_id;
    vector<double> gen_rate; // 每顆衛星對該需求的生成率
    Requirement() : id(0), gs1(0), gs2(0), served_by_sat_id(-1) {}
    Requirement(int _id, int _gs1, int _gs2, int numSats)
        : id(_id), gs1(_gs1), gs2(_gs2), served_by_sat_id(-1), gen_rate(numSats, 0.0) {}
};

// 給優先佇列使用的結構，依生成率排序（越大優先）
struct Requirement_queue {
    int gs1, gs2, sat;
    int req_id;
    double gen_rate;
    Requirement_queue(int _gs1, int _gs2, int _sat, int _req_id, double _gen_rate)
        : gs1(_gs1), gs2(_gs2), sat(_sat), req_id(_req_id), gen_rate(_gen_rate) {}
};

struct CompareRequirement {
    bool operator()(const Requirement_queue &a, const Requirement_queue &b) const {
        return a.gen_rate < b.gen_rate;
    }
};

Satellite sat[1000];
Requirement req[1000];

// 資料前處理，建立每顆衛星可服務的地面站集合以及計算每個需求在各衛星下的生成率
void data_process(){
    // 對每個地面站、每顆衛星檢查 fidelity 是否達標
    for (int i = 0; i < G; i++){
        for (int j = 0; j < S; j++){
            if (sat[j].check_fid(i)){
                sat[j].gs_serve.insert(i);
            }
        }
    }
    
    // 對每個需求，計算每顆衛星能否同時服務該需求的兩個地面站，
    // 若可以則該需求在此衛星下的生成率為兩地面站生成率的較小值
    for (int i = 0; i < R; i++){
        req[i].gen_rate.resize(S, 0.0);
        for (int j = 0; j < S; j++){
            if (sat[j].check_both_gs(req[i].gs1, req[i].gs2)){
                sat[j].gsp_serve.insert(i);
                double rate1 = rate_gs_sat[{req[i].gs1, j}];
                double rate2 = rate_gs_sat[{req[i].gs2, j}];
                req[i].gen_rate[j] = min(rate1, rate2);
            }
        }
    }
}

// 貪婪演算法：利用優先佇列從高生成率開始分配，限制每顆衛星只能服務一個需求
void greedy_algorithm(){
    priority_queue<Requirement_queue, vector<Requirement_queue>, CompareRequirement> req_queue;
    for (int i = 0; i < R; i++){
        for (int j = 0; j < S; j++){
            if (req[i].gen_rate[j] == 0.0) continue;
            req_queue.push(Requirement_queue(req[i].gs1, req[i].gs2, j, i, req[i].gen_rate[j]));
        }
    }
    set<int> sat_available;
    for (int i = 0; i < S; i++){
        sat_available.insert(i);
    }
    
    while (!req_queue.empty() && !sat_available.empty()){
        auto cur = req_queue.top();
        req_queue.pop();
        if (sat_available.find(cur.sat) == sat_available.end() || req[cur.req_id].served_by_sat_id != -1){
            continue;
        }
        req[cur.req_id].served_by_sat_id = cur.sat;
        sat_available.erase(cur.sat);
    }
}

// 只讀取必要的資料：需求與每個衛星與地面站間的 fidelity 與 generation rate
// 輸入檔 dataset.txt 格式如下：
// 第一行: S G R
// 接下來 R 行: 每行兩個數字 (gs1 gs2)
// 接著依照每顆衛星、每個地面站的順序，讀取兩個數字：fidelity 與 generation rate
void input(){
    ifstream in("dataset.txt");
    assert(in);
    in >> S >> G >> R;
    // 讀取需求：兩個地面站編號
    for (int i = 0; i < R; i++){
        in >> req[i].gs1 >> req[i].gs2;
        req[i].id = i;
        req[i].gen_rate.resize(S, 0.0);
    }
    // 讀取每顆衛星與每個地面站間的資料（僅兩個參數）
    for (int i = 0; i < S; i++){
        for (int j = 0; j < G; j++){
            double tmp_fid, tmp_gen_rate;
            in >> tmp_fid >> tmp_gen_rate;
            fid_gs_sat[{j, i}] = tmp_fid;
            rate_gs_sat[{j, i}] = tmp_gen_rate;
        }
    }
}

// 輸出總生成率：累計所有已分配需求在其服務衛星下的生成率
void output(){
    ofstream out("res_greedy.txt");
    assert(out);
    double total = 0;
    for (int i = 0; i < R; i++){
        int sat_id = req[i].served_by_sat_id;
        if (sat_id != -1){
            total += req[i].gen_rate[sat_id];
            out << "req " << i << " is served by " << sat_id << " with gen rate " << req[i].gen_rate[sat_id] << '\n';
        }
    }
    out << "Total generation rate: " << total << "\n";
    out.close();
}

int main(){
    input();
    data_process();
    greedy_algorithm();
    output();
    return 0;
}
