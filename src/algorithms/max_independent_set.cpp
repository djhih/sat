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

Satellite sat[20000];
GroundStation gs[20000];
Requirement req[200000];

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
                // if(req[i].gs1 == 0 && req[i].gs2 == 2 && j == 2){
                // cout << "rate1 " << rate1 << " rate2 " << rate2 << '\n';
                // cout << rate_gs_p_sat[{0, 2}] << '\n';}
                req[i].gen_rate[j] = (min(rate1, rate2));
            }
        }
    }
}

struct Node {
    int gs1, gs2, gsp_id, sat, id;
    bool erased = 0;
    double weight_degree;
    vector<int>neighboor;
    Node() : gs1(-1), gs2(-1), gsp_id(-1), sat(-1), weight_degree(-1), id(-1), neighboor(vector<int>()) {}
    Node(int gs1, int gs2, int gsp_id, int sat) : 
        gs1(gs1), gs2(gs2), gsp_id(gsp_id), sat(sat), weight_degree(-1), neighboor(vector<int>()) {}
};

vector<Node>nodes;

void print_nodes(){
    for(int i=0; i<nodes.size(); i++){
        cout << "node " << i << '\n';
        cout << "gs1 " << nodes[i].gs1 << " gs2 " << nodes[i].gs2 << " sat " << nodes[i].sat << '\n';
        cout << "weight degree " << nodes[i].weight_degree << '\n';
        cout << "neighboor: ";
        for(auto x:nodes[i].neighboor) cout << x << ' ';
        cout << '\n';
    }
}

void transfer_graph(){
    int ptr = 0;
    for (int i=0; i<S; i++) {
        auto gss = sat[i].gsp_serve;
        for(auto x:gss){
            int tmp_gs1 = req[x].gs1, tmp_gs2 = req[x].gs2;
            if(tmp_gs1 > tmp_gs2)
                swap(tmp_gs1, tmp_gs2); 
            Node tmp(tmp_gs1, tmp_gs2, x, i); // gs1, gs2, req_id, sat
            tmp.id = ptr;
            nodes.push_back(tmp); 
            ptr++;
        }
    }  
    cout << "<MIS>: finish transfer_graph and get nodes size " << nodes.size() << '\n';
    for(int i=0; i<ptr; i++){
        for(int j=i+1; j<ptr; j++){
            if(nodes[i].gs1 == nodes[j].gs1 || nodes[i].gs2 == nodes[j].gs2 || nodes[i].gs1 == nodes[j].gs2 || nodes[i].gs2 == nodes[j].gs1 || nodes[i].sat == nodes[j].sat){
                nodes[i].neighboor.emplace_back(j);
                nodes[j].neighboor.emplace_back(i);
            }
        }
    }
    for(int i=0; i<ptr; i++){
        double weight = 0;
        for(auto v: nodes[i].neighboor){
            weight += rate_gs_p_sat[{nodes[v].gsp_id, nodes[v].sat}];
        }
        weight /= (double)(nodes[i].neighboor.size());
        nodes[i].weight_degree = weight;
    }
    // print_nodes();
}

vector<int>ans;
void dfs(){
    // 1. find minimum weight_degree node to add to answer
    // 2. erase the node's neighboor
    int can_use_count = nodes.size();
    assert(nodes.size() < 1000000);
    while(can_use_count > 0){
        int cur_choose = -1, cur_wd = 10000000;
        for(int i=0; i<nodes.size(); i++){
            if(nodes[i].erased) continue;
            if(nodes[i].weight_degree < cur_wd){
                cur_wd = nodes[i].weight_degree;
                cur_choose = i;
            }
        }
        if(cur_choose != -1)
            ans.push_back(cur_choose);
        // cout << "push " << cur_choose << '\n';
        // cut
        nodes[cur_choose].erased = 1;
        can_use_count--;
        // cout << "erase node choose " << cur_choose << '\n';
        for(auto nid:nodes[cur_choose].neighboor){
            nodes[nid].erased = 1;
            // cout << "erase node " << nid << '\n';
            can_use_count--;
        }
    }
    
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


string outfile = "dataset/output/res_max.txt";
void output(){
    ofstream out(outfile);
    assert(out);
    double tot_rate = 0;
    //accept gs1 1 gs2 2 sat 1 gen rate 27
    for(int i=0; i<ans.size(); i++){
        out << "accept gs1 " << nodes[ans[i]].gs1 << " gs2 " << nodes[ans[i]].gs2 
            << " sat " << nodes[ans[i]].sat << " gen rate " << rate_gs_p_sat[{nodes[ans[i]].gsp_id, nodes[ans[i]].sat}] << '\n';
        tot_rate += rate_gs_p_sat[{nodes[ans[i]].gsp_id, nodes[ans[i]].sat}];
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
    cout << "finish input\n";
    data_process();
    cout << "finish data process\n";
    transfer_graph();
    cout << "finish transfer graph\n";
    dfs();
    cout << "finish dfs\n";
    output();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Total elapsed time: " << elapsed.count() << " seconds" << std::endl;
    return 0;
}
