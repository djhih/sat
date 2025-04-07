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
#include <filesystem>
#include <direct.h>
#include <limits.h>
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
Requirement req[1000000];

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

void greedy_algorithm(){
    priority_queue<Requirement_queue, vector<Requirement_queue>, CompareRequirement> req_queue;
    for (int i = 0; i < R; i++){
        for (int j = 0; j < S; j++){
            if (req[i].gen_rate[j] == -1) continue;
            Requirement_queue tmp(req[i].gs1, req[i].gs2, j, i, req[i].gen_rate[j]);
            req_queue.push(tmp); 
            // cout << "push " << req[i].gen_rate[j] << '\n';
        }
    }
    set<int> sat_can_use, gs_can_use;
    for (int i = 0; i < S; i++){
        sat_can_use.insert(i);
    }
    for (int i=0; i<G; i++){
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
        sat_can_use.erase(tp.sat);
        gs_can_use.erase(tp.gs1);
        gs_can_use.erase(tp.gs2);
    }
}

string infile = "dataset/raw/dataset.txt";
void input(){
    ifstream in(infile);
    char buffer[1000];
    if (_getcwd(buffer, 1000)) {
        std::cout << "Current working directory: " << buffer << std::endl;
    } else {
        perror("_getcwd error");
    }
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


string outfile = "dataset/output/res_greedy.txt";
void output(){
    ofstream out(outfile);
    assert(out);
    double ans = 0;
    //accept gs1 1 gs2 2 sat 1 gen rate 27
    for (int i = 0; i < R; i++){
        int cur_served_sat_id = req[i].served_by_sat_id;
        if (cur_served_sat_id != -1 && req[i].gen_rate[cur_served_sat_id] != 0){
            ans += req[i].gen_rate[cur_served_sat_id];
            out << "accept gs1 " << req[i].gs1 << " gs2 " << req[i].gs2 << " sat " << cur_served_sat_id << " gen rate " << req[i].gen_rate[cur_served_sat_id] << '\n';
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
    greedy_algorithm();
    output();
    return 0;
}
