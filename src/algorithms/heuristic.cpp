#define _USE_MATH_DEFINES 
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <queue>
#include <map>
#include <set>
#include <fstream>
#include <cassert>
#include <unordered_map>
using namespace std;

typedef pair<double, double> pdd;
typedef pair<int, int> pii;

// 自訂 hash 結構
struct PairHash {
    size_t operator()(const pii& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};

// 全域變數，使用 unordered_map 搭配自訂 hash
int S, G, R;
unordered_map<pii, double, PairHash> dis_gs_sat;       // dis(gs_i, sat_j)
unordered_map<pii, double, PairHash> ang_gs_sat;       // elevation angle
unordered_map<pii, double, PairHash> h_atm_gs_sat;     // distance from ground to atmosphere (gs-sat line)
unordered_map<pii, bool, PairHash> served_gs_sat;      
unordered_map<pii, bool, PairHash> served_gs_p_sat;    
unordered_map<pii, double, PairHash> rate_gs_sat;      
unordered_map<pii, double, PairHash> rate_gs_p_sat;    
unordered_map<pii, double, PairHash> fid_gs_sat;       
unordered_map<pii, double, PairHash> fid_gs_p_sat;       

double fidelity_threshold = 0.5;
double angle_threshold = 20;
double alpha;

struct ECEF {
    double x, y, z;
    ECEF() : x(0.0), y(0.0), z(0.0) {}
    ECEF(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
};

struct Satellite {
    int id;
    ECEF location;
    set<int> gs_serve;  // 可以服務的 ground station 集合
    set<int> gsp_serve; // 可以服務的 ground station pair 集合
    int cur_serve_req_id;
    
    Satellite() : id(-1), location(), cur_serve_req_id(-1), gs_serve(set<int>()), gsp_serve(set<int>()) {}  
    Satellite(const int id, const ECEF &loc, const set<int> &gs_set, const set<int> &gsp_set, int served)
        : id(id), location(loc), gs_serve(gs_set), gsp_serve(gsp_set), cur_serve_req_id(served) {}
    
    bool check_ang(int gs) {
        return ang_gs_sat[{gs, id}] <= angle_threshold;
    }
    
    bool check_fid(int gs) {
        return fid_gs_sat[{gs, id}] >=  fidelity_threshold;
    }
    
    bool check_both_gs(int gs1, int gs2) {
        return (check_fid(gs1)) && (check_fid(gs2));
    }
};

struct GroundStation {
    int id;
    ECEF location;
    
    GroundStation() : id(0), location() {}
    GroundStation(int _id, const ECEF &loc) : id(_id), location(loc) {}
};

struct Requirement {
    int id;
    int gs1, gs2;
    int served_by_sat_id;
    vector<double> gen_rate; // generation rate for each satellite

    Requirement() : id(0), gs1(0), gs2(0), served_by_sat_id(-1), gen_rate(vector<double>()) {}
    Requirement(int _id, int _gs1, int _gs2, int numSats)
        : id(_id), gs1(_gs1), gs2(_gs2), served_by_sat_id(-1), gen_rate(numSats, 0.0) {}
};

struct Requirement_queue {
    int gs1, gs2, sat;
    int req_id;
    double gen_rate;
    
    Requirement_queue() : gs1(-1), gs2(-1), sat(-1), req_id(-1), gen_rate(-1) {}
    Requirement_queue(int _gs1, int _gs2, int _sat, int _req_id, double _gen_rate)
        : gs1(_gs1), gs2(_gs2), sat(_sat), req_id(_req_id), gen_rate(_gen_rate) {}
};

struct CompareRequirement {
    bool operator()(const Requirement_queue &a, const Requirement_queue &b) const {
        return a.gen_rate < b.gen_rate;
    }
};

Satellite sat[10000];
GroundStation gs[10000];
Requirement req[100000];

void data_process(){
    
    for (int i = 0; i < G; i++){
        for (int j = 0; j < S; j++){
            if (sat[j].check_fid(i)){
                sat[j].gs_serve.emplace(i);
            }
        }
    }
    
    for (int i = 0; i < R; i++){
        req[i].gen_rate.resize(S);
        for (int j = 0; j < S; j++){
            if (sat[j].check_both_gs(req[i].gs1, req[i].gs2)){
                sat[j].gsp_serve.emplace(i);
                double rate1 = rate_gs_sat[{req[i].gs1, j}];
                double rate2 = rate_gs_sat[{req[i].gs2, j}];
                rate_gs_p_sat[{i, j}] = min(rate1, rate2);
                // cout << "rate1 " << rate1 << " rate2 " << rate2 << '\n';
                req[i].gen_rate[j] = (min(rate1, rate2));
            }
        }
    }
}

bool dfs(int reqq, vector<bool>& visited, vector<int>& match, const vector<vector<int>>& graph) {
    for (int sat_id : graph[reqq]) {
        if (!visited[sat_id]) {
            visited[sat_id] = true;
            if (match[sat_id] == -1 || dfs(match[sat_id], visited, match, graph)) {
                match[sat_id] = reqq;
                return true;
            }
        }
    }
    return false;
}
// 用貪婪法求最大加權匹配（近似最佳解）
// 此版本將所有 (requirement, satellite) 組合（來自 rate_gs_p_sat）
// 根據 generation rate 由大到小排序，再依序挑選符合條件的配對。
void efficient_weighted_matching() {
    // 建立所有候選的配對：元素格式 (gen_rate, requirement id, satellite id)
    vector<tuple<double, int, int>> edges;
    for (auto &entry : rate_gs_p_sat) {
        int req_id = entry.first.first;   // requirement 編號
        int sat_id = entry.first.second;    // satellite 編號
        double gen_rate = entry.second;
        // 只考慮 generation rate 大於 0 的情況
        if (gen_rate > 0)
            edges.push_back(make_tuple(gen_rate, req_id, sat_id));
    }
    
    // 依 generation rate 由大到小排序
    sort(edges.begin(), edges.end(), [](const auto &a, const auto &b) {
        return get<0>(a) > get<0>(b);
    });
    
    // 建立記錄各個 requirement、satellite 與 ground station 是否已被使用的標記
    vector<bool> req_used(R, false);
    vector<bool> sat_used(S, false);
    vector<bool> gs_used(G, false);
    
    // 逐一考慮最高權重的邊
    for (auto &edge : edges) {
        double rate;
        int req_id, sat_id;
        tie(rate, req_id, sat_id) = edge;
        // 若此 requirement 已被配對則跳過
        if (req_used[req_id])
            continue;
        
        // 取得此 requirement 牽涉的兩個 ground station
        int gs1 = req[req_id].gs1;
        int gs2 = req[req_id].gs2;
        // 若該 requirement 的兩個 ground station與該衛星都還未被使用，
        // 則進行配對：設定該 requirement 由此 sat 服務，
        // 並標記 requirement、satellite 與兩個 ground station均已被使用
        if (!gs_used[gs1] && !gs_used[gs2] && !sat_used[sat_id]) {
            req[req_id].served_by_sat_id = sat_id;
            // 同時記錄該 candidate 的 generation rate（可用於後續統計）
            req[req_id].gen_rate[sat_id] = rate;
            req_used[req_id] = true;
            sat_used[sat_id] = true;
            gs_used[gs1] = true;
            gs_used[gs2] = true;
        }
    }
    
    // 對未成功配對的 requirement，統一設為 -1
    for (int i = 0; i < R; i++) {
        if (!req_used[i])
            req[i].served_by_sat_id = -1;
    }
}

void Hungarian_algo(){
    // 左側為 requirement，右側為 satellite
    int n_rows = R; // requirement 數量
    int n_cols = S; // satellite 數量
    int n = max(n_rows, n_cols); // 建立 n x n 的矩陣
    // 建立成本矩陣：cost[i][j] = - generation rate (若有邊)
    // 若該對 (i, j) 沒有對應 edge (或 generation rate 為 0)，則設為 0
    vector<vector<double>> cost(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if(i < n_rows && j < n_cols){
                pii key = {i, j};
                double rate = 0.0;
                if(rate_gs_p_sat.find(key) != rate_gs_p_sat.end()){
                    rate = rate_gs_p_sat[key];
                }
                cost[i][j] = -rate; // 轉換成成本，因為我們要最大化 generation rate
            }
            else{
                cost[i][j] = 0;
            }
        }
    }
    
    // 以下為標準的 Hungarian 演算法 (1-index版本)
    const double INF = 1e9;
    vector<double> u(n+1, 0), v(n+1, 0);
    vector<int> p(n+1, 0), way(n+1, 0);
    for (int i = 1; i <= n; i++){
        p[0] = i;
        int j0 = 0;
        vector<double> minv(n+1, INF);
        vector<bool> used(n+1, false);
        do {
            used[j0] = true;
            int i0 = p[j0];
            int j1 = 0;
            double delta = INF;
            for (int j = 1; j <= n; j++){
                if(!used[j]){
                    double cur = cost[i0-1][j-1] - u[i0] - v[j];
                    if(cur < minv[j]){
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if(minv[j] < delta){
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            for (int j = 0; j <= n; j++){
                if(used[j]){
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while(p[j0] != 0);
        do{
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while(j0);
    }
    
    vector<int> row_match(n_rows, -1);
    for (int j = 1; j <= n; j++){
        int i = p[j];
        if(i - 1 < n_rows && j - 1 < n_cols){
            row_match[i-1] = j - 1;
        }
    }
    
    for (int i = 0; i < n_rows; i++){
        int sat_id = row_match[i];
        pii key = {i, sat_id};
        double rate = 0.0;
        if(sat_id != -1 && rate_gs_p_sat.find(key) != rate_gs_p_sat.end()){
            rate = rate_gs_p_sat[key];
        }
        if(sat_id != -1 && rate > 0){
            req[i].served_by_sat_id = sat_id;
        }
        else {
            req[i].served_by_sat_id = -1;
        }
    }
}

void greedy_algorithm(){
    // 建立所有可能的 candidate：每個元素為 (gen_rate, requirement id, satellite id)
    vector<tuple<double, int, int>> assignments;
    for (auto &entry : rate_gs_p_sat) {
        int req_id = entry.first.first;   // requirement pair 編號
        int sat_id = entry.first.second;    // satellite 編號
        double gen_rate = entry.second;
        // 只考慮 generation rate > 0 的配對
        if (gen_rate > 0) {
            assignments.push_back(make_tuple(gen_rate, req_id, sat_id));
        }
    }
    
    // 依 generation rate 由大到小排序
    sort(assignments.begin(), assignments.end(), [](const auto &a, const auto &b) {
        return get<0>(a) > get<0>(b);
    });
    
    // 追蹤每個 requirement、satellite 以及 ground station 是否已配對
    vector<bool> req_assigned(R, false);
    vector<bool> sat_assigned(S, false);
    vector<bool> gs_assigned(G, false);
    double total_gen_rate = 0.0;
    
    // 依序選取最高權重的 candidate
    for (auto &assignment : assignments) {
        double gen_rate;
        int req_id, sat_id;
        tie(gen_rate, req_id, sat_id) = assignment;
        
        // 若該 requirement 已配對，則跳過
        if (req_assigned[req_id])
            continue;
        
        int gs1 = req[req_id].gs1;
        int gs2 = req[req_id].gs2;
        // 若兩個 ground station 和該衛星都還沒被使用，則進行配對
        if (!gs_assigned[gs1] && !gs_assigned[gs2] && !sat_assigned[sat_id]) {
            req[req_id].served_by_sat_id = sat_id;
            // 記錄此配對的 generation rate（覆蓋原本 req[i].gen_rate[sat_id] 也可）
            req[req_id].gen_rate[sat_id] = gen_rate;
            // cout << "req " << req_id << " sat " << sat_id << " gen rate " << gen_rate << '\n';
            gs_assigned[gs1] = true;
            gs_assigned[gs2] = true;
            sat_assigned[sat_id] = true;
            req_assigned[req_id] = true;
            total_gen_rate += gen_rate;
        }
    }
    
    // 對於未配對的 requirement，設為 -1
    for (int i = 0; i < R; i++){
        if (!req_assigned[i])
            req[i].served_by_sat_id = -1;
    }
    
    // cout << "Total generation rate by greedy: " << total_gen_rate << endl;
}

string infile = "dataset/raw/dataset.txt";
void input(){
    ifstream in(infile);
    assert(in);
    in >> S >> G >> R;
    // cout << "S " << S << " G " << G << " R " << R << '\n';
    for (int i = 0; i < S; i++){
        in >> sat[i].location.x >> sat[i].location.y >> sat[i].location.z;
        sat[i].id = i;
        // cout << sat[i].location.x << ' ' << sat[i].location.y << ' ' << sat[i].location.z << '\n';
    }
    for (int i = 0; i < G; i++){
        in >> gs[i].location.x >> gs[i].location.y >> gs[i].location.z;
        gs[i].id = i;
    }
    for (int i = 0; i < R; i++){
        int a, b;
        in >> a >> b;
        if(a>b) swap(a, b);
        req[i].gs1 = a;
        req[i].gs2 = b;
        req[i].id = i;
        req[i].gen_rate.resize(S, 0.0);
    }
    for (int i = 0; i < S; i++){
        for (int j = 0; j < G; j++){
            int tmp_dis, tmp_ang, tmp_h_atm;
            double tmp_fid, tmp_gen_rate;
            in >> tmp_dis >> tmp_ang >> tmp_h_atm >> tmp_fid >> tmp_gen_rate;
            dis_gs_sat[{j, i}] = tmp_dis;
            ang_gs_sat[{j, i}] = tmp_ang;
            h_atm_gs_sat[{j, i}] = tmp_h_atm;
            fid_gs_sat[{j, i}] = tmp_fid;
            rate_gs_sat[{j, i}] = tmp_gen_rate;
        }
    }
}


string outfile = "dataset/output/res_he.txt";
void output(){
    ofstream out(outfile);
    assert(out);
    double ans = 0;
    //accept gs1 1 gs2 2 sat 1 gen rate 27
    for (int i = 0; i < R; i++){
        int cur_served_sat_id = req[i].served_by_sat_id;
        if (cur_served_sat_id != -1 && req[i].gen_rate[cur_served_sat_id] != 0){
            ans += req[i].gen_rate[cur_served_sat_id];
            out << "accept " << "gs1 " << req[i].gs1 << " gs2 " << req[i].gs2 << " sat " << cur_served_sat_id << " gen rate " << req[i].gen_rate[cur_served_sat_id] << '\n';
        }
    }
    out << "Total generation rate: " << ans << '\n';
    out.close();
}

int main(int argc, char* argv[]){
    if(argc == 3){
        infile = argv[1];
        outfile = argv[2];
    }
    input();
    data_process();
    // Hungarian_algo();
    efficient_weighted_matching();
    greedy_algorithm();
    output();
    return 0;
}
