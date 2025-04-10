#pragma O3
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
#include <omp.h>
using namespace std;

typedef pair<double, double> pdd;
typedef pair<int, int> pii;

// 自訂 hash 結構
struct PairHash {
    size_t operator()(const pii& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};

// 全域變數與資料結構
int S, G, R;
unordered_map<pii, double, PairHash> dis_gs_sat;    // dis(gs_i, sat_j)
unordered_map<pii, double, PairHash> ang_gs_sat;      // elevation angle
unordered_map<pii, double, PairHash> h_atm_gs_sat;    // ground-to-atmosphere distance
unordered_map<pii, bool, PairHash> served_gs_sat;     
unordered_map<pii, bool, PairHash> served_gs_p_sat;   
unordered_map<pii, double, PairHash> rate_gs_sat;     
unordered_map<pii, double, PairHash> rate_gs_p_sat;   
unordered_map<pii, double, PairHash> fid_gs_sat;      
unordered_map<pii, double, PairHash> fid_gs_p_sat;    

// 使用固定大小的二維陣列記錄節點間是否有連結（若節點數量過多需調整）
bool node_adj[1000][1000];    
// 使用固定大小的布林陣列記錄某些節點是否已在改善解中
bool nodes_in_imp_ans[1000000]; 

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
    set<int> gs_serve;  // 可服務的 ground station 集合
    set<int> gsp_serve; // 可服務的 ground station pair 集合
    int cur_serve_req_id;
    
    Satellite() : id(-1), location(), cur_serve_req_id(-1),
                  gs_serve(set<int>()), gsp_serve(set<int>()) {}  
    Satellite(const int id, const ECEF &loc, const set<int> &gs_set, const set<int> &gsp_set, int served)
        : id(id), location(loc), gs_serve(gs_set), gsp_serve(gsp_set), cur_serve_req_id(served) {}
    
    bool check_ang(int gs) {
        return ang_gs_sat[{gs, id}] <= angle_threshold;
    }
    
    bool check_fid(int gs) {
        return fid_gs_sat[{gs, id}] >= fidelity_threshold;
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
    int served_by_sat_id; // 用於貪婪求解
    vector<double> gen_rate; // 各衛星產生速率

    Requirement() : id(0), gs1(0), gs2(0), served_by_sat_id(-1),
                    gen_rate(vector<double>()) {}
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

Satellite sat[10009];
GroundStation gs[10900];
Requirement req[100000];

// 進行資料處理，計算每個衛星可以服務哪些 ground station / pair
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
                req[i].gen_rate[j] = min(rate1, rate2);
            }
        }
    }
}

struct Node {
    int gs1, gs2, gsp_id, sat, id;
    bool erased = 0;
    double weight;
    vector<int> neighboor;
    Node() : gs1(-1), gs2(-1), gsp_id(-1), sat(-1), weight(-1), id(-1),
             neighboor(vector<int>()) {}
    Node(int gs1, int gs2, int gsp_id, int sat) : 
        gs1(gs1), gs2(gs2), gsp_id(gsp_id), sat(sat), weight(-1),
        neighboor(vector<int>()) {}
    void print(){
        cout << "gs1 " << gs1 << " gs2 " << gs2 << " sat " << sat << " weight " << weight << "\n";
    }
};

vector<Node> nodes;

// 建立圖，將每個衛星之可服務 ground station pair 轉換成圖中的節點，並依據條件建立鄰接關係
void transfer_graph(){
    int ptr = 0;
    for (int i = 0; i < S; i++) {
        auto gss = sat[i].gsp_serve;
        for(auto x : gss){
            int tmp_gs1 = req[x].gs1, tmp_gs2 = req[x].gs2;
            if(tmp_gs1 > tmp_gs2)
                swap(tmp_gs1, tmp_gs2); 
            Node tmp(tmp_gs1, tmp_gs2, x, i); // gs1, gs2, req_id, sat
            tmp.id = ptr;
            nodes.push_back(tmp); 
            ptr++;
        }
    }  
    for(int i = 0; i < ptr; i++){
        for(int j = i+1; j < ptr; j++){
            if(nodes[i].gs1 == nodes[j].gs1 || nodes[i].gs2 == nodes[j].gs2 ||
               nodes[i].gs1 == nodes[j].gs2 || nodes[i].gs2 == nodes[j].gs1 ||
               nodes[i].sat == nodes[j].sat){
                nodes[i].neighboor.push_back(j);
                nodes[j].neighboor.push_back(i);
                node_adj[i][j] = 1;
            }
        }
    }
    for(int i = 0; i < ptr; i++){
        sort(nodes[i].neighboor.begin(), nodes[i].neighboor.end());
    }
    for(int i = 0; i < ptr; i++){
        nodes[i].weight = rate_gs_p_sat[{nodes[i].gsp_id, nodes[i].sat}]; 
    }

    ofstream out("dataset/output/graph.txt");
    assert(out);
    out << "nodes size " << nodes.size() << '\n';
    for(int i = 0; i < nodes.size(); i++){
        out << "node " << i << '\n';
        out << "gs1 " << nodes[i].gs1 << " gs2 " << nodes[i].gs2 
            << " sat " << nodes[i].sat << '\n';
        out << "weight " << nodes[i].weight << '\n';
        out << "neighboor ";
        for(auto x : nodes[i].neighboor){
            out << x << ' ';
        }
        out << '\n';
    }
    out.close();
}

vector<pair<int, int>> greedy_ans_gsp_sat;
void greedy(){
    priority_queue<Requirement_queue, vector<Requirement_queue>, CompareRequirement> req_queue;
    for (int i = 0; i < R; i++){
        for (int j = 0; j < S; j++){
            if (req[i].gen_rate[j] == -1) continue;
            Requirement_queue tmp(req[i].gs1, req[i].gs2, j, i, req[i].gen_rate[j]);
            req_queue.push(tmp); 
        }
    }
    set<int> sat_can_use, gs_can_use;
    for (int i = 0; i < S; i++){
        sat_can_use.insert(i);
    }
    for (int i = 0; i < G; i++){
        gs_can_use.insert(i);
    }
    while (!req_queue.empty() && !sat_can_use.empty()){
        auto tp = req_queue.top(); 
        req_queue.pop();
        if ((sat_can_use.count(tp.sat) == 0 || gs_can_use.count(tp.gs1) == 0 || gs_can_use.count(tp.gs2) == 0)
            || req[tp.req_id].served_by_sat_id != -1){
            continue;
        }
        int cur_rq_id = tp.req_id;
        req[cur_rq_id].served_by_sat_id = tp.sat;
        greedy_ans_gsp_sat.push_back({cur_rq_id, tp.sat});
        sat_can_use.erase(tp.sat);
        gs_can_use.erase(tp.gs1);
        gs_can_use.erase(tp.gs2);
    }
}

set<int> imp_ans_nodes;

void rescale(double val_greedy){
    const double K = 200, V = nodes.size();
    double ratio = K * V / val_greedy;
    for(int i = 0; i < V; i++){
        // cout << "nodes " << i << " weight " << nodes[i].weight << " * " << ratio << " = " << '\n';
        nodes[i].weight *= ratio;
    }
}

// 函式 check_improve 用來檢查傳入的 claws 是否能改善解答
bool check_improve(vector<vector<int>>& claws){
    // 對每個 claw 收集刪除候選節點 (del_nodes)
    vector<vector<int>> del_nodes;
    for(auto claw : claws){
        set<int> del_tmp;
        for(auto x : claw){
            for(auto y : nodes[x].neighboor){
                if(nodes_in_imp_ans[y]){
                    del_tmp.insert(y);
                }
            }
        }
        del_nodes.push_back(vector<int>(del_tmp.begin(), del_tmp.end()));
    }
    // 依序檢查每個 claw 的增益是否大於門檻值
    bool found = false;
    for(int _claw = 0; _claw < claws.size(); _claw++){
        auto claw = claws[_claw];
        vector<int> del_tmp = del_nodes[_claw];
        double add_val = 0;
        for(auto x : claw){
            if(nodes_in_imp_ans[x]){
                cout << "already in ans " << x << '\n';
                continue;
            }
            add_val += nodes[x].weight * nodes[x].weight;
        }
        for(auto x : del_tmp){
            add_val -= nodes[x].weight * nodes[x].weight;
        }
        // cout << "add_val " << add_val << '\n';
        if(add_val > 1e-6){
            // 若改善成立，加入新解且移除受影響節點
            for(auto x : claw){ 
                imp_ans_nodes.insert(x);
                nodes_in_imp_ans[x] = 1;
            }
            for(auto x : del_tmp){
                imp_ans_nodes.erase(x);
                nodes_in_imp_ans[x] = 0;
            }
            cout << "add claw " << '\n';
            for(auto x : claw){
                cout << "gs1 " << nodes[x].gs1 << " gs2 " << nodes[x].gs2 
                     << " sat " << nodes[x].sat << '\n';
            }
            cout << "del claw " << '\n';
            for(auto x : del_tmp){
                cout << "gs1 " << nodes[x].gs1 << " gs2 " << nodes[x].gs2 
                     << " sat " << nodes[x].sat << '\n';
            }
            cout << "add_val " << add_val << '\n';
            found = true;
            return found;
        }
    }
    return found;
}

// 平行化 get_claw 函式：針對中心節點，在其鄰居中平行搜尋各種候選 claw 組合
bool get_claw(int center) {
    int nei_cen = nodes[center].neighboor.size();
    vector<int> reduce_nodes;
    for (auto x : nodes[center].neighboor) {
        if (!nodes_in_imp_ans[x])
            reduce_nodes.push_back(x);
    }
    int reduce_size = reduce_nodes.size();
    bool found = false;
    
    // 平行化三個節點組合搜尋
    #pragma omp parallel for schedule(dynamic) shared(found)
    for (int i0 = 0; i0 < reduce_size; i0++) {
        if (found) continue;
        int i = reduce_nodes[i0];
        if (nodes_in_imp_ans[i]) continue;
        for (int j0 = i0 + 1; j0 < reduce_size; j0++) {
            if (found) break;
            int j = reduce_nodes[j0];
            if (nodes_in_imp_ans[j]) continue;
            for (int k0 = j0 + 1; k0 < reduce_size; k0++) {
                if (found) break;
                int k = reduce_nodes[k0];
                if (nodes_in_imp_ans[k]) continue;
                vector<vector<int>> claws;
                if (!node_adj[i][j] && !node_adj[i][k] && !node_adj[j][k]) {
                    claws.push_back({i, j, k});
                }
                if (!claws.empty() && check_improve(claws)) {
                    #pragma omp critical
                    {
                        found = true;
                    }
                    break;
                }
            }
        }
    }
    if (found)
        return true;
    
    // 平行化兩個節點的候選 (claw2)
    vector<vector<int>> claw2;
    #pragma omp parallel for schedule(dynamic)
    for (int i0 = 0; i0 < reduce_size; i0++) {
        int i = reduce_nodes[i0];
        if (nodes_in_imp_ans[i]) continue;
        for (int j0 = i0 + 1; j0 < reduce_size; j0++) {
            int j = reduce_nodes[j0];
            if (nodes_in_imp_ans[j] || node_adj[i][j]) continue;
            #pragma omp critical
            {
                claw2.push_back({i, j});
            }
        }
    }
    if (check_improve(claw2))
        return true;
    
    // 平行化單一節點候選 (claw1)
    vector<vector<int>> claw1;
    #pragma omp parallel for schedule(dynamic)
    for (int i0 = 0; i0 < reduce_size; i0++) {
        int i = reduce_nodes[i0];
        if (nodes_in_imp_ans[i]) continue;
        #pragma omp critical
        {
            claw1.push_back({i});
        }
    }
    return check_improve(claw1);
}

void imp(){
    // Step 1: 透過貪婪法取得初始解
    greedy();
    double greedy_value = 0;
    for(auto x : greedy_ans_gsp_sat){
        auto [gsp, sat] = x;
        greedy_value += rate_gs_p_sat[{gsp, sat}];
    }
    rescale(greedy_value);
    // 將貪婪初始解複製到 imp_ans_nodes 與 nodes_in_imp_ans 陣列中
    for(int i = 0; i < greedy_ans_gsp_sat.size(); i++){
        for(int j = 0; j < nodes.size(); j++){
            if(nodes[j].gsp_id == greedy_ans_gsp_sat[i].first && nodes[j].sat == greedy_ans_gsp_sat[i].second){
                imp_ans_nodes.insert(j);
            }
        }
    }
    for(auto x : imp_ans_nodes){
        nodes_in_imp_ans[x] = 1;
    }
    // Step 2: 針對每個節點當作中心搜尋改善解 (claw)
    for(int i = 0; i < nodes.size(); i++){
        // cout << "now on " << i << '\n';
        int center = i;
        if(get_claw(center)){
            // 若發現改善解則重置 i 以重新搜尋
            cout << "found claw " << center << '\n';
            i = -1;
        }
    }    
}

string infile = "dataset/raw/dataset4.txt";
void input(){
    ifstream in(infile);
    assert(in);
    in >> S >> G >> R;
    for (int i = 0; i < S; i++){
        in >> sat[i].location.x >> sat[i].location.y >> sat[i].location.z;
        sat[i].id = i;
    }
    for (int i = 0; i < G; i++){
        in >> gs[i].location.x >> gs[i].location.y >> gs[i].location.z;
        gs[i].id = i;
    }
    for (int i = 0; i < R; i++){
        int a, b;
        in >> a >> b;
        if(a > b) swap(a, b);
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
    cout << "end input\n";
}
string outfile = "dataset/output/res_imp_4.txt";
void output(){
    ofstream out(outfile);
    assert(out);
    double tot_rate = 0;
    for(auto x : imp_ans_nodes){
        out << "accept gs1 " << nodes[x].gs1 << " gs2 " << nodes[x].gs2 
            << " sat " << nodes[x].sat << " gen rate " 
            << rate_gs_p_sat[{nodes[x].gsp_id, nodes[x].sat}] << '\n';
        tot_rate += rate_gs_p_sat[{nodes[x].gsp_id, nodes[x].sat}];
    }
    out << "Total generation rate: " << tot_rate << '\n';
    out.close();
}

int main(int argc, char* argv[]){
    if(argc == 3){
        infile = argv[1];
        outfile = argv[2];
    }
    input();
    data_process();
    transfer_graph();
    imp();
    output();
    return 0;
}
