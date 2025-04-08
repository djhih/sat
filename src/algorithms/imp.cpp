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
unordered_map<pii, bool, PairHash> node_adj;    

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
    int served_by_sat_id; //! use for greedy
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

Satellite sat[10009];
GroundStation gs[10900];
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
    double weight;
    vector<int>neighboor;
    Node() : gs1(-1), gs2(-1), gsp_id(-1), sat(-1), weight(-1), id(-1), neighboor(vector<int>()) {}
    Node(int gs1, int gs2, int gsp_id, int sat) : 
        gs1(gs1), gs2(gs2), gsp_id(gsp_id), sat(sat), weight(-1), neighboor(vector<int>()) {}

    void print(){
        cout << "gs1 " << gs1 << " gs2 " << gs2 << " sat " << sat << " weight " << weight << "\n";
    }
};

vector<Node>nodes;

// void print_nodes(){
//     for(int i=0; i<nodes.size(); i++){
//         cout << "node " << i << '\n';
//         cout << "gs1 " << nodes[i].gs1 << " gs2 " << nodes[i].gs2 << " sat " << nodes[i].sat << '\n';
//         cout << "weight degree " << nodes[i].weight << '\n';
//         cout << "neighboor: ";
//         for(auto x:nodes[i].neighboor) cout << x << ' ';
//         cout << '\n';
//     }
// }

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
    for(int i=0; i<ptr; i++){
        for(int j=i+1; j<ptr; j++){
            if(nodes[i].gs1 == nodes[j].gs1 || nodes[i].gs2 == nodes[j].gs2 || nodes[i].gs1 == nodes[j].gs2 || nodes[i].gs2 == nodes[j].gs1 || nodes[i].sat == nodes[j].sat){
                nodes[i].neighboor.emplace_back(j);
                nodes[j].neighboor.emplace_back(i);
                node_adj[{i, j}] = 1;
            }
        }
    }
    for(int i=0; i<ptr; i++){
        sort(nodes[i].neighboor.begin(), nodes[i].neighboor.end());
    }

    for(int i=0; i<ptr; i++){
        nodes[i].weight = rate_gs_p_sat[{nodes[i].gsp_id, nodes[i].sat}]; 
    }
    //print_nodes();
}

vector<pair<int, int>>greedy_ans_gsp_sat;
void greedy(){
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
        req[cur_rq_id].served_by_sat_id = tp.sat;  // not really needed
        greedy_ans_gsp_sat.push_back({cur_rq_id, tp.sat});
        sat_can_use.erase(tp.sat);
        gs_can_use.erase(tp.gs1);
        gs_can_use.erase(tp.gs2);
    }
}

// void print_greedy_ans(){
//     cout << "greedy ans " << '\n';
//     for(auto x: greedy_ans_gsp_sat){
//         auto [gsp, sat] = x;
//         cout << "gs1 " << req[gsp].gs1 << " gs2 " << req[gsp].gs2 << " sat " << sat << '\n';
//     }
// }
set<int>imp_ans_nodes;

void rescale(double val_greedy){
    const double K = 200, V = nodes.size();
    double ratio = K*V/val_greedy;

    for(int i=0; i<V; i++){
        nodes[i].weight *= ratio;
    }
}

// DFS to find all claws
// We need to count add_val in DFS
// Biggest claw has 3 talons(=3)
// we will set cur_claw as the talons we plan to add, if we find 3 talons, then we check if it can improve the answer
// if it can improve, then we will add all the talons into the answer, remove the node connected to the talons, and return with value 1
// if it cannot improve, then we will return with value 0

set<int> last_claw;
bool DFS(int center, int cnt, vector<bool> &visited, set<int>cur_claw){
    if(cnt > 0){
        // chekc if it can improve the answer
        // 1. get all delete nodes
        // 2. check if it can improve the answer
        // 3. if it can improve, then add all, and remove the node connected to the talons

        set<int>del_tmp;
        for(auto x:cur_claw){
            for(auto y:nodes[x].neighboor){
                if(imp_ans_nodes.count(y) == 1){
                    del_tmp.emplace(y);
                }
            }
        }

        double add_val = 0;
        for(auto x: cur_claw){
            add_val += nodes[x].weight * nodes[x].weight;
        }
        for(auto x: del_tmp){
            add_val -= nodes[x].weight * nodes[x].weight;
        }
        
        if(add_val > 0){
            // add claw
            for(auto x: cur_claw){ 
                imp_ans_nodes.emplace(x);
            }
            // remove del_tmp
            for(auto x: del_tmp){
                imp_ans_nodes.erase(x);
            }

            cout << "add claw " << '\n';
            for(auto x: cur_claw){
                cout << "gs1 " << nodes[x].gs1 << " gs2 " << nodes[x].gs2 << " sat " << nodes[x].sat << '\n';
            }
            cout << "del claw " << '\n';
            for(auto x: del_tmp){
                cout << "gs1 " << nodes[x].gs1 << " gs2 " << nodes[x].gs2 << " sat " << nodes[x].sat << '\n';
            }
            cout << "add_val " << add_val << '\n';
            return 1;
        }

        return 0;
    }

    for(int i=0; i<nodes[center].neighboor.size(); i++){
        int nei = nodes[center].neighboor[i];
        if(visited[nei] == 0){
            visited[nei] = 1;

            bool connected = 0;
            // check if nei connect to the cur_claw
            for(auto x: nodes[nei].neighboor){
                if(cur_claw.count(x)){
                    connected = 1;
                    break;
                }
            }
            if(connected == 1){
                continue;
            }
            
            cur_claw.emplace(nei);
            DFS(center, cnt+1, visited, cur_claw);
            cur_claw.erase(nei);
            DFS(center, cnt+1, visited, cur_claw);
        }
    }
    return 0;
}

// loop version
bool get_claw(int center){
    int nei_cen = nodes[center].neighboor.size();
    for(int i=0; i<nei_cen; i++){
        if(imp_ans_nodes.count(i) == 1){
            continue;
        }
        for(int j=i+1; j<nei_cen; j++){
            if(imp_ans_nodes.count(j) == 1){
                continue;
            }
            for(int k=j+1; k<nei_cen; k++){
                if(imp_ans_nodes.count(k) == 1){
                    continue;
                }
                // check if the three nodes are connected to the node in answer set
                // check if thses talons are connected with each other, 
                //      if so, split to multiple claws, and store into vector
                // cout << "i " << i << " j " << j << " k " << k << '\n';
                vector<vector<int>> claws;
                bool connected12 = 0, connected13 = 0, connected23 = 0;
                if(node_adj.count({i, j}) == 1){
                    connected12 = 1;
                }
                if(node_adj.count({i, k}) == 1){
                    connected13 = 1;
                }
                if(node_adj.count({j, k}) == 1){
                    connected23 = 1;
                }
                
                if(!connected12 && !connected13 && !connected23){
                    claws.push_back({i, j, k});
                }
                if(!connected12){
                    claws.push_back({i, j});
                }
                if(!connected13){
                    claws.push_back({i, k});
                }
                if(!connected23){
                    claws.push_back({j, k});
                }
                

                // find del nodes
                vector<vector<int>> del_nodes;
                for(auto claw: claws){
                    vector<int>del_tmp;
            
                    for(auto x: claw){
                        for(auto y: nodes[x].neighboor){
                            if(imp_ans_nodes.count(y) == 1){
                                del_tmp.emplace_back(y);
                            }
                        }
                    }
                    del_nodes.push_back(del_tmp);
                }

                // check if the claw can improve the answer
                bool found = 0;
                for(int _claw = 0; _claw < claws.size(); _claw++){
                    auto claw = claws[_claw];
    
                    vector<int>del_tmp = del_nodes[_claw];
                    double add_val = 0;
                    for(auto x: claw){
                        if(imp_ans_nodes.count(x) == 1){
                            cout << "already in ans " << x << '\n';
                            continue;
                        }
                        add_val += nodes[x].weight * nodes[x].weight;
                    }
                    for(auto x: del_tmp){
                        add_val -= nodes[x].weight * nodes[x].weight;
                    }
                    
                    if(add_val > 0){
                        // add claw
                        for(auto x: claw){ 
                            imp_ans_nodes.emplace(x);
                        }
                        // remove del_tmp
                        for(auto x: del_tmp){
                            imp_ans_nodes.erase(x);
                        }

                        cout << "add claw " << '\n';
                        for(auto x: claw){
                            cout << "gs1 " << nodes[x].gs1 << " gs2 " << nodes[x].gs2 << " sat " << nodes[x].sat << '\n';
                        }
                        cout << "del claw " << '\n';
                        for(auto x: del_tmp){
                            cout << "gs1 " << nodes[x].gs1 << " gs2 " << nodes[x].gs2 << " sat " << nodes[x].sat << '\n';
                        }
                        cout << "add_val " << add_val << '\n';
                        found = 1;
                        break;
                    } else {
                        cout << "not add claw " << ' ';
                        cout << "add_val " << add_val << ' ';
                        cout << "i " << i << " j " << j << " k " << k << ' ';
                        cout << "gs1 " << nodes[claw[0]].gs1 << " gs2 " << nodes[claw[0]].gs2 << " sat " << nodes[claw[0]].sat << '\n';
                    }
                }
                if(found == 1){
                    return found;
                }
               
            }
        }
    }

    vector<vector<int>> claws;
    for(int i=0; i<nei_cen; i++){
        if(imp_ans_nodes.count(i) == 1){
            continue;
        }
        claws.push_back({i});
    }
    vector<vector<int>> del_nodes;
    for(auto claw: claws){
        vector<int>del_tmp;
        for(auto x: claw){
            for(auto y: nodes[x].neighboor){
                if(imp_ans_nodes.count(y) == 1){
                    del_tmp.emplace_back(y);
                }
            }
        }
        del_nodes.push_back(del_tmp);
    }
    bool found = 0;
    for(int _claw = 0; _claw < claws.size(); _claw++){
        auto claw = claws[_claw];
        vector<int>del_tmp = del_nodes[_claw];
        double add_val = 0;
        for(auto x: claw){
            add_val += nodes[x].weight * nodes[x].weight;
        }
        for(auto x: del_tmp){
            add_val -= nodes[x].weight * nodes[x].weight;
        }
        if(add_val > 0){
            // add claw
            for(auto x: claw){ 
                imp_ans_nodes.emplace(x);
            }
            // remove del_tmp
            for(auto x: del_tmp){
                imp_ans_nodes.erase(x);
            }

            cout << "add claw " << '\n';
            for(auto x: claw){
                cout << "gs1 " << nodes[x].gs1 << " gs2 " << nodes[x].gs2 << " sat " << nodes[x].sat << '\n';
            }
            cout << "del claw " << '\n';
            for(auto x: del_tmp){
                cout << "gs1 " << nodes[x].gs1 << " gs2 " << nodes[x].gs2 << " sat " << nodes[x].sat << '\n';
            }
            cout << "add_val " << add_val << '\n';
            found = 1;
            break;
        } else {
            cout << "not add claw " << '\n';
        }
    }
    if(found == 1){
        return found;
    }
    return 0;
}


void imp(){
    // 1. greedy to get inital answer
    // 2. find all claws in graph(nodes)
    // 3. everytime we find a claw, we have to check if it can improve answer
    //    3.1 how to check?
    //    3.2 we have to removed the node connected to the talons
    // 4. if it can improve, then add all, and remove the node connected to the talons  
    //    4.1 after add ans, goto 2.
    // 5. if not go to next claw, and do 3.2 until no more claw can add
    
    greedy(); // 1. get greedy_ans_gsp_sat
    // print_greedy_ans();

    double greedy_value = 0;
    for(auto x: greedy_ans_gsp_sat){
        auto [gsp, sat] = x;
        greedy_value += rate_gs_p_sat[{gsp, sat}];
    }

    rescale(greedy_value);

    // print_greedy_ans();

    // copy answer
    for(int i=0; i<greedy_ans_gsp_sat.size(); i++){
        for(int j=0; j<nodes.size(); j++){
            if(nodes[j].gsp_id == greedy_ans_gsp_sat[i].first && nodes[j].sat == greedy_ans_gsp_sat[i].second){
                imp_ans_nodes.emplace(j);
            }
        }
    }

    // 2. find all claws
    for(int i=0; i<nodes.size(); i++){
        cout << "now on " << i << '\n';
        int center = i, add_val = 0;     // choose nodes[i] as center node
        bool loop_res = get_claw(center);
        if(loop_res == 1){
            i = -1;
        }
    }    
    
}

string infile = "dataset/raw/dataset1.txt";
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
    cout << "end input\n";
}
string outfile = "dataset/output/res_imp_1.txt";
void output(){
    ofstream out(outfile);
    assert(out);
    double tot_rate = 0;
    for(auto x:imp_ans_nodes){
        out << "accept gs1 " << nodes[x].gs1 << " gs2 " << nodes[x].gs2 
            << " sat " << nodes[x].sat << " gen rate " << rate_gs_p_sat[{nodes[x].gsp_id, nodes[x].sat}] << '\n';
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
    // dfs();
    imp();
    output();
    return 0;
}
