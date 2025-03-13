#define _USE_MATH_DEFINES 
#include<algorithm>
#include<iostream>
#include<vector>
#include<cmath>
#include<queue>
#include<map>
#include<set>
using namespace std;

typedef pair<double, double> pdd;
typedef pair<int, int> pii;

struct Hash {
    size_t operator()(const pii& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};
struct Hash2 {
    size_t operator()(const pii& p) const {
        size_t h1 = hash<int>()(p.first);
        size_t h2 = hash<int>()(p.second);
        return h1 ^ (h2 << 1); // 使用 XOR 和位移減少碰撞
    }
};

int S, G, R;

// global variable soooooooooo baddddddddddddd 
map<pii, double, Hash>dis_gs_sat;       // dis_gs_sat[{i, j}] = dis(gs_i, sat_j)
map<pii, double, Hash>ang_gs_sat;       // = the elervation of gs_i and sat_j
map<pii, double, Hash>h_atm_gs_sat;     // = distance from ground to atmosphere (gs-sat line)
map<pii, bool, Hash>served_gs_sat;      // = if sat_j can serve gs_i;
map<pii, bool, Hash>served_gs_p_sat;    // = if sat_j can serve gsp_i
map<pii, double, Hash>rate_gs_sat;      // generation rate for gs_i and sat_j
map<pii, double, Hash>rate_gs_p_sat;    // generation rate for gsp_i and sat_j
map<pii, double, Hash>fid_gs_sat;       // fidelity for gs_i and sta_j
map<pii, double, Hash>fid_gs_p_sat;     // fidelity for gs_i and sat_j

double fidelity_threshold = 0.5;
double angle_threshold = 20;
double alpha;

struct ECEF {
    double x, y, z;
    ECEF() : x(0.0), y(0.0), z(0.0) {}
    ECEF(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

};

/* not in use, we set fidelity in dataset generator
double get_fid(double Fi0, double nj, double nij) {
    return 0.25 * (1 + (4 * Fi0 - 1) / pow(1 + nj / nij, 2));
}

double get_eta(double dT_i, double dR_k, double lambda_ik, double alpha, double h_ik) {
    double pi = M_PI;  
    double term1 = pi * pow(dT_i / 2, 2);
    double term2 = pi * pow(dR_k / 2, 2);
    double denominator = pow(lambda_ik, 2);
    double exp_term = exp(-alpha * h_ik);
    return (term1 * term2 / denominator) * exp_term;
}
*/

struct Satellite {
    int id;
    ECEF location;
    set<int>gs_serve; // the gs this sat can serve
    set<int>gsp_serve; // the gs pair this sat can serve
    int cur_serve_req_id;

    Satellite() : id(-1), location(), cur_serve_req_id(-1), gs_serve(set<int>()), gsp_serve(set<int>())  {}  
    Satellite(const int id, const ECEF &loc, const std::set<int> &gs_set, const std::set<int> &gsp_set, int served)
        : id(id), location(loc), gs_serve(gs_set), gsp_serve(gsp_set), cur_serve_req_id(served) {}

    bool check_ang(int gs){
        return ang_gs_sat[{gs, id}] <= angle_threshold;
    }

    bool check_fid(int req){
        return fid_gs_p_sat[{req, id}] >= fidelity_threshold;
    }

    bool check_both_gs(int gs1, int gs2){
        bool res1 = (check_ang(gs1) & check_fid(gs1));
        bool res2 = (check_ang(gs2) & check_fid(gs2));
        return res1 & res2;
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
    vector<double> gen_rate; // generation rate  //! not init, don't forget constructor

    Requirement() : id(0), gs1(), gs2(), served_by_sat_id(-1), gen_rate(R, -1) {}
    Requirement(int _id, const int _gs1, const int _gs2, int numSats)
        : id(_id), gs1(_gs1), gs2(_gs2), served_by_sat_id(-1), gen_rate(numSats, 0.0) {}

    void set_gen_rate(int sat_id){
        //! if sat cannot serve this requirement, set the val to -1
        if(rate_gs_sat[{gs1, sat_id}] == 0 || rate_gs_sat[{gs2, sat_id}] == 0){
            gen_rate[sat_id] = -1;
            rate_gs_p_sat[{id, sat_id}] = -1;
        } else {
            gen_rate[sat_id] = min(rate_gs_sat[{gs1, sat_id}], rate_gs_sat[{gs2, sat_id}]);
            rate_gs_p_sat[{id, sat_id}] = gen_rate[sat_id];  
        }
    }
};

struct Requirement_queue{ // all id
    int gs1, gs2, sat;
    int req_id;
    double gen_rate;

    Requirement_queue() : gs1(-1), gs2(-1), sat(-1), req_id(-1), gen_rate(-1) {}
    Requirement_queue(int _gs1, int _gs2, int _sat, int _req_id, double _gen_rate)
        : gs1(_gs1), gs2(_gs2), sat(_sat), req_id(_req_id), gen_rate(_gen_rate) {}
    
};

struct CompareRequirement {
    bool operator()(const Requirement_queue &a, const Requirement_queue &b) const {
        return a.gen_rate < b.gen_rate; // gen_rate 較高的將會排在前面
    }
};

Satellite sat[100];
GroundStation gs[100];
Requirement req[100];

void greedy_algorithm(){
    // everytime choose the best rate req to served
    // 1. put all req into queue 
    // 2. assign the maximum gen_rate to sat
    // 3. if already accept or all sat has been used then stop

    priority_queue<Requirement_queue>req_queue;
    // init req queue
    for(int i=0; i<R; i++){
        for(int j=0; j<S; j++){
            if(req[i].gen_rate[j] == -1) continue;
            Requirement_queue tmp(req[i].gs1, req[i].gs2, j, i, req[i].gen_rate[j]);
            req_queue.push(tmp); 
        }
    }
    // init sat_can_use set
    set<bool>sat_can_use;
    for(int i=0; i<S; i++){
        sat_can_use.insert(i);
    }

    while(!req_queue.empty() && !sat_can_use.empty()){
        auto tp = req_queue.top(); 
        req_queue.pop();
        if(sat_can_use.count(tp.sat) == 0 || req[tp.req_id].served_by_sat_id != -1){
            continue;
        }
        int cur_rq_id = tp.req_id;
        req[cur_rq_id].served_by_sat_id = tp.sat; 
        sat_can_use.erase(tp.sat);
    }
}

void input(){
    cin >> S >> G >> R;
    for(int i=0; i<S; i++){
        cin >> sat[i].location.x >> sat[i].location.y >> sat[i].location.z;
        sat[i].id = i;
        int tmp_gs, tmp_dis, tmp_ang, tmp_h_atm, tmp_fid, tmp_gen_rate;
        cin >> tmp_gs >> tmp_dis >> tmp_ang >> tmp_h_atm >> tmp_fid >> tmp_gen_rate;
        dis_gs_sat[{tmp_gs, i}] = tmp_dis;
        ang_gs_sat[{tmp_gs, i}] = tmp_ang;
        h_atm_gs_sat[{tmp_gs, i}] = tmp_h_atm;
        fid_gs_sat[{tmp_gs, i}] = tmp_fid;
        rate_gs_sat[{tmp_gs, i}] = tmp_gen_rate;
    }
    for(int i=0; i<G; i++){
        cin >> gs[i].location.x >> gs[i].location.y >> gs[i].location.z;
        gs[i].id = i;
    }
    for(int i=0l; i<R; i++){
        cin >> req[i].gs1 >> req[i].gs2;
        req[i].id = i;
    }
}

void data_proccess(){
    // add gs to sat served list
    for(int i=0; i<G; i++){
        for(int j=0; j<S; j++){
            if(sat[j].check_fid(i) && sat[j].check_ang(i)){
                sat[j].gs_serve.emplace(i);
            }
        }
    }

    // add gs pait to sat served list
    for(int i=0; i<R; i++){
        for(int j=0; j<S; j++){
            if(sat[j].check_both_gs(req[i].gs1, req[i].gs2)){
                sat[j].gsp_serve.emplace(i);

                double rate1 = rate_gs_sat[{req[i].gs1, j}];
                double rate2 = rate_gs_sat[{req[i].gs2, j}];
                rate_gs_p_sat[{i, j}] = min(rate1, rate2);
            }
        }
    }
}

void output(){
    // add all ac req gen rate
    double ans = 0;
    for(int i=0; i<R; i++){
        int cur_served_sat_id = req[i].served_by_sat_id;
        if(cur_served_sat_id != -1){
            ans += req[i].gen_rate[cur_served_sat_id];
        }
    }
    cout << "tot rate: " << ans << '\n';
}

int main(){
    // input
    // data_proccess
    // greedy_algorithm
    // end
}