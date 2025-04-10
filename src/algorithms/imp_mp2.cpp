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
#include <chrono>
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

////////////////////////
// 舊程式碼的部分 (data_process, transfer_graph, greedy, input, output) 保持不變
////////////////////////

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
        nodes[i].weight *= ratio;
    }
}

// *********************************************************************
// 以下新增「批次搜尋候選改善」的機制
// *********************************************************************

// 結構 Candidate：記錄候選改善包含新增的節點、要刪除的節點，以及淨增益值
struct Candidate {
    vector<int> add_nodes;
    vector<int> remove_nodes;
    double gain;
};

// evaluate_candidate() 計算傳入 claw 的淨增益
Candidate evaluate_candidate(const vector<int>& claw) {
    Candidate cand;
    cand.add_nodes = claw;
    cand.gain = 0.0;
    set<int> removals;
    for (auto x : claw) {
        // 若 x 未在解中，計算加成
        if (!nodes_in_imp_ans[x])
            cand.gain += nodes[x].weight * nodes[x].weight;
        // 將 x 的所有鄰居中在解中的節點加入 removals 集合
        for (auto y : nodes[x].neighboor) {
            if (nodes_in_imp_ans[y])
                removals.insert(y);
        }
    }
    cand.remove_nodes.assign(removals.begin(), removals.end());
    double removalPenalty = 0.0;
    for (auto y : cand.remove_nodes) {
        removalPenalty += nodes[y].weight * nodes[y].weight;
    }
    cand.gain -= removalPenalty;
    return cand;
}

// 批次搜尋所有中心節點的候選改善
bool batch_candidate_improvement() {
    vector<Candidate> candidates;  // 收集所有正增益候選
    // 平行遍歷每個 center node
    #pragma omp parallel for schedule(dynamic)
    for (int center = 0; center < (int)nodes.size(); center++) {
        vector<int> reduce_nodes;
        cout << "now on " << center << '\n';
        // 對 center node，取得其鄰居中尚未在解中的節點
        for (auto x : nodes[center].neighboor) {
            if (!nodes_in_imp_ans[x])
                reduce_nodes.push_back(x);
        }
        int reduce_size = reduce_nodes.size();
        // 枚舉三個節點組合（claw3）
        for (int i = 0; i < reduce_size; i++){
            for (int j = i+1; j < reduce_size; j++){
                for (int k = j+1; k < reduce_size; k++){
                    int a = reduce_nodes[i], b = reduce_nodes[j], c = reduce_nodes[k];
                    if (!node_adj[a][b] && !node_adj[a][c] && !node_adj[b][c]) {
                        vector<int> claw = {a, b, c};
                        Candidate cand = evaluate_candidate(claw);
                        if (cand.gain > 1e-6) {
                            #pragma omp critical
                            candidates.push_back(cand);
                        }
                    }
                }
            }
        }
        // 枚舉兩個節點組合（claw2）
        for (int i = 0; i < reduce_size; i++){
            for (int j = i+1; j < reduce_size; j++){
                int a = reduce_nodes[i], b = reduce_nodes[j];
                if (!node_adj[a][b]) {
                    vector<int> claw = {a, b};
                    Candidate cand = evaluate_candidate(claw);
                    if (cand.gain > 1e-6) {
                        #pragma omp critical
                        candidates.push_back(cand);
                    }
                }
            }
        }
        // 單一節點候選（claw1）
        for (int i = 0; i < reduce_size; i++){
            vector<int> claw = {reduce_nodes[i]};
            Candidate cand = evaluate_candidate(claw);
            if (cand.gain > 1e-6) {
                #pragma omp critical
                candidates.push_back(cand);
            }
        }
    } // end parallel for

    // 若沒有候選改善，則回傳 false
    if(candidates.empty())
        return false;

    // 選擇淨增益最大的候選改善（你也可以選第一個正增益候選）
    Candidate bestCand = candidates[0];
    for(auto &cand : candidates) {
        if(cand.gain > bestCand.gain)
            bestCand = cand;
    }

    // 在單一執行緒中更新全域解：將 bestCand.add_nodes 加入解，並將 bestCand.remove_nodes 移除
    for(auto x : bestCand.add_nodes){
        imp_ans_nodes.insert(x);
        nodes_in_imp_ans[x] = true;
    }
    for(auto x : bestCand.remove_nodes){
        imp_ans_nodes.erase(x);
        nodes_in_imp_ans[x] = false;
    }

    cout << "Applied candidate: gain = " << bestCand.gain << endl;
    return true;
}

// 修改 imp() 函數，採用批次更新策略
void imp() {
    // Step 1: 透過貪婪法取得初始解
    greedy();
    double greedy_value = 0;
    for(auto x : greedy_ans_gsp_sat){
        auto [gsp, sat] = x;
        greedy_value += rate_gs_p_sat[{gsp, sat}];
    }
    rescale(greedy_value);
    // 將貪婪初始解複製到 imp_ans_nodes 與 nodes_in_imp_ans 陣列中
    for (int i = 0; i < greedy_ans_gsp_sat.size(); i++){
        for (int j = 0; j < nodes.size(); j++){
            if(nodes[j].gsp_id == greedy_ans_gsp_sat[i].first && nodes[j].sat == greedy_ans_gsp_sat[i].second){
                imp_ans_nodes.insert(j);
            }
        }
    }
    for(auto x : imp_ans_nodes){
        nodes_in_imp_ans[x] = true;
    }
    // Step 2: 批次搜尋候選改善，直到找不到任何改善為止
    while(batch_candidate_improvement()){
        // 每次更新後可以印出目前解的狀況
        cout << "Batch improvement applied, current solution size: " << imp_ans_nodes.size() << '\n';
    }
}

// 以下 input() 與 output() 與 main() 保持原有內容
string infile = "dataset/raw/dataset1.txt";
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
string outfile = "dataset/output/res_imp_mp_1.txt";
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
    auto start = std::chrono::high_resolution_clock::now();
    data_process();
    transfer_graph();
    imp();
    output();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Total elapsed time: " << elapsed.count() << " seconds" << std::endl;
    return 0;
}
