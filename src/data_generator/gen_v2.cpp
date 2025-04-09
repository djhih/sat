#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <functional>
#include <utility>
#include <set>
#include <cassert>
#include <random>
#include <time.h>
#include <filesystem>
#include <limits.h>
#include <string>
using namespace std;

const double PI = acos(-1.0);
int num_sat = 20;
int num_gs = 30;              
int num_req = 900; 

// --- vector operation ---
double dot(const vector<double>& a, const vector<double>& b) {
    double dp = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        dp += a[i] * b[i];
    }
    return dp;
}

vector<double> subtract(const vector<double>& a, const vector<double>& b) {
    return { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
}

vector<double> scalar_divide(const vector<double>& v, double s) {
    return { v[0]/s, v[1]/s, v[2]/s };
}

double norm(const std::vector<double>& v) {
    return std::sqrt(dot(v, v));
}

std::vector<double> normalize(const std::vector<double>& v) {
    double n = norm(v);
    if(n == 0) throw std::runtime_error("Zero vector cannot be normalized");
    std::vector<double> u(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        u[i] = v[i] / n;
    }
    return u;
}

// phy & coordinate conversion
// (R, theta, phi) -> (x,y,z)
vector<double> generate_point(double R, double theta, double phi) { 
    return { R * sin(phi) * cos(theta),
             R * sin(phi) * sin(theta),
             R * cos(phi) };
}

double get_R_earth() {
    return 6.371e6; 
}

// --- ��z�v�����p�� ---
double eta_sg(double L, double lam = 810e-9, double w0 = 0.025, double rg = 0.75) {
    double LR = PI * w0 * w0 / lam;
    double wg = w0 * sqrt(1 + (L*L) / (LR*LR));
    return 1 - exp(-2 * rg * rg / (wg * wg));
}

double eta_Atm(double L, double h, double etaZen = 0.5) {
    double R = get_R_earth();
    double denominator = (h / L) - ((L * L - h * h) / (2 * R * L));
    double secZen = 1.0 / denominator;
    return pow(etaZen, secZen);
}

double eta_Tot(double L, double h) {
    return eta_sg(L) * eta_Atm(L, h);
}

// --- fidelity ---
double fidelity_Fij(double theta_k1, double theta_e, double n_sg, double n_ij, double F0) {
    ofstream out("debug_fid.txt", std::ios::app);

    if ((theta_k1) >= theta_e) {
        out << "nsg " << n_sg << " n_ij " << n_ij << " F0 " << F0 << '\n';
        out << "Fid " << 0.25 * (1.0 + (4.0 * F0 - 1.0) / pow((1.0 + n_sg / n_ij), 2)) << '\n';
        out.close();
        return 0.25 * (1.0 + (4.0 * F0 - 1.0) / pow((1.0 + n_sg / n_ij), 2));
    } 
        
    out.close();
    return 0.0;
}

//! not sure how to decide gen rate 
double gen_rate(double d) {
    return 1;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis2(0.01, 1);
    // uniform_int_distribution<> dis2(1, 10);
    return dis2(gen);
}

// --- �s���Z���P�̨ΰ��� ---
double link_distance(double d, double h) {
    double R = get_R_earth();
    return sqrt(4.0 * R * (R + h) * pow(sin(d / (4.0 * R)), 2) + h * h);
}

// --- ���� ---
double compute_elevation_angle(const vector<double>& ground_ecef, const vector<double>& sat_ecef, bool degrees = false) {
    vector<double> diff = subtract(sat_ecef, ground_ecef);
    double norm_diff = norm(diff);
    vector<double> up = scalar_divide(ground_ecef, norm(ground_ecef));
    double dot_up = dot(diff, up);
    double elev_rad = asin(dot_up / norm_diff);
    if (degrees)
        return elev_rad * 180.0 / PI;
    else
        return elev_rad;
}

struct GS_Sat_Data {
    int dis;
    int ang;
    int h_atm;
    double fid;
    double gen_rate;
};

vector<vector<double>> satellites;
vector<vector<double>> ground_stations;

double ground_to_atmosphere_distance(const std::vector<double>& gs_ecef, const std::vector<double>& sat_ecef) {
    double R = 6371000;    
    double atmosphere_thickness = 10000;
    std::vector<double> d(3);
    for (size_t i = 0; i < 3; i++) {
        d[i] = sat_ecef[i] - gs_ecef[i];
    }
    d = normalize(d);  // ���o���V�q

    double R_atm = R + atmosphere_thickness;

    double gs_dot_d = dot(gs_ecef, d);
    double gs_norm_sq = dot(gs_ecef, gs_ecef);
    double A = 1.0;
    double B = 2.0 * gs_dot_d;
    double C = gs_norm_sq - R_atm * R_atm;

    double discriminant = B * B - 4 * A * C;
    if (discriminant < 0) {
        cerr << "gs ecef " << gs_ecef[0] << " " << gs_ecef[1] << " " << gs_ecef[2] << '\n';
        cerr << "sat ecef " << sat_ecef[0] << " " << sat_ecef[1] << " " << sat_ecef[2] << '\n';
        cerr << "gs norm " << gs_norm_sq << '\n';
        throw std::runtime_error("No intersection: ���Ӥ��|�o�͡A�]���a�������b�j��h��");
    }

    double sqrt_disc = std::sqrt(discriminant);
    // �D�o��Ӯ� t1 �M t2�A�q�` t1 ���t�]�I�V�g�u��V�^�At2 ����
    double t1 = (-B - sqrt_disc) / (2 * A);
    double t2 = (-B + sqrt_disc) / (2 * A);

    // ��� t �����B�̤p���ȡ]�� gs �̪񪺥��I�^
    double t = (t1 >= 0) ? t1 : t2;
    if (t < 0) {
        throw std::runtime_error("�g�u�P�j��h�L���V���I");
    }

    return t; // t ���q�a�����u�g�u�ܥ��I���Z�� (m)
}


string sat_filename = "dataset/code/output/satellite_coordinates_selected.txt";
void input_sat_data(vector<vector<double>>& satellites) {
    ifstream in(sat_filename);
    assert(in);
    in >> num_sat;;
    satellites.resize(num_sat, vector<double>(3));
    for (int i = 0; i < num_sat; i++) {
        in >> satellites[i][0] >> satellites[i][1] >> satellites[i][2];
    }
}

string gs_filename = "dataset/code/output/gs_loc.txt";
void input_GS_data(vector<vector<double>>& ret) {
    ifstream in(gs_filename);
    char buffer[1000];
    assert(in);  
    set<vector<double>> gs_sat;
    in >> num_gs;
    for (int i = 0; i < num_gs; i++) {
        vector<double> tmp(3);  
        in >> tmp[0] >> tmp[1] >> tmp[2];
        gs_sat.insert(tmp);
    }
    ret = vector<vector<double>>(gs_sat.begin(), gs_sat.end());
}

vector<pair<int, int>> random_ground_satation_pair(int num, int G_total) {
    vector<pair<int, int>> result;
    set<pair<int, int>> res;
    
    if (G_total < 2) {
        cerr << "Error: G_total must be at least 2!\n";
        return result; // �Ū� vector
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, num_gs - 1);

    set<pair<int, int>>pair_tried;
    while(res.size() < num && pair_tried.size() < G_total * G_total){
        int a = dis(gen);
        int b = dis(gen);
        if(a > b) swap(a, b);
        while(a == b || res.count({b, a}) || pair_tried.count({a, b})){
            a = dis(gen);
            b = dis(gen);
            if(a > b) swap(a, b);
        }
        pair_tried.insert({a, b});
        // check if there exist a sat can serve gs a and gs b, which means fidelity of (gsa, sati) and (gsb, sati) > 0.5
        bool can_serve = 0;
        int idx = 0;
        for(auto sat: satellites){
            double tmp_dis_a = norm(subtract(sat, ground_stations[a]));
            double tmp_dis_b = norm(subtract(sat, ground_stations[b]));
            double tmp_ang_a = compute_elevation_angle(ground_stations[a], sat, true);
            double tmp_ang_b = compute_elevation_angle(ground_stations[b], sat, true);
            double tmp_h_atm_a = ground_to_atmosphere_distance(ground_stations[a], sat);
            assert(tmp_h_atm_a > 0);
            double tmp_h_atm_b = ground_to_atmosphere_distance(ground_stations[b], sat);
            assert(tmp_h_atm_b > 0);
            double fid0 = 0.85;  
            double height = norm(sat) - get_R_earth();
            double tmp_fid_a = fidelity_Fij(tmp_ang_a, 20.0, eta_Tot(tmp_dis_a, tmp_h_atm_a), 1.0, fid0);
            double tmp_fid_b = fidelity_Fij(tmp_ang_b, 20.0, eta_Tot(tmp_dis_b, tmp_h_atm_b), 1.0, fid0);
            if(tmp_fid_a > 0.5 && tmp_fid_b > 0.5){
                can_serve = 1;
                break;
            }
            idx++;
        }
        if(can_serve){
            res.insert({a, b});
            // cout << a << " and " << b << " caqn served by " << idx << '\n';
        } 
            
    }

    result = vector<pair<int, int>>(res.begin(), res.end());
    cerr << "try " << pair_tried.size() << " pairs\n";
    cerr << "get " << result.size() << " pairs\n";
    return result; 
}

vector<pair<int, int>> get_all_ground_station_pairs(int lim, int GS) {
    vector<pair<int, int>> pairs;
    int cnt = 0;
    for (int i = 0; i < GS && cnt < lim; i++) {
        for (int j = i + 1; j < GS && cnt < lim; j++) {
            bool can_serve = 0;
            for(auto sat: satellites){
                double tmp_dis_a = norm(subtract(sat, ground_stations[i]));
                double tmp_dis_b = norm(subtract(sat, ground_stations[j]));
                double tmp_ang_a = compute_elevation_angle(ground_stations[i], sat, true);
                double tmp_ang_b = compute_elevation_angle(ground_stations[j], sat, true);
                double tmp_h_atm_a = ground_to_atmosphere_distance(ground_stations[i], sat);
                assert(tmp_h_atm_a > 0);
                double tmp_h_atm_b = ground_to_atmosphere_distance(ground_stations[j], sat);
                assert(tmp_h_atm_b > 0);
                double fid0 = 0.85;  
                double height = norm(sat) - get_R_earth();
                double tmp_fid_a = fidelity_Fij(tmp_ang_a, 20.0, eta_Tot(tmp_dis_a, tmp_h_atm_a), 1.0, fid0);
                double tmp_fid_b = fidelity_Fij(tmp_ang_b, 20.0, eta_Tot(tmp_dis_b, tmp_h_atm_b), 1.0, fid0);
                if(tmp_fid_a > 0.6 && tmp_fid_b > 0.6){
                    can_serve = 1;
                    break;
                }
            }
            if(can_serve){
                pairs.push_back({i, j});
                cnt++;
            }
                
        }
    }
    return pairs;
}


int main(int argc, char* argv[]) {
    ofstream out("debug_fid.txt");
    out.clear();
    out.close();
      
    double R_earth_val = get_R_earth();   // �a�y�b�| (m)
    double h = 500e3;             // �ìP���� (m)
    
    string output_filename = "dataset/raw/dataset.txt";

    if(argc == 5){
        num_req = stoi(argv[1]);
        sat_filename = argv[2];
        gs_filename = argv[3];
        output_filename = argv[4];
    } 

    // --- Sat loc ---
    input_sat_data(satellites);
    cerr << satellites.size() << " get all sat data\n";

    // --- GS loc ---
    input_GS_data(ground_stations);
    cerr << "get all gs data\n";

    // --- GSP ---
    vector<pair<int, int>> requirements = get_all_ground_station_pairs(num_req, num_gs);
    // vector<pair<int, int>> requirements = random_ground_satation_pair(num_req, num_gs);
    cerr << "get all gs pairs\n";

    // vector<pair<int, int>> requirements = generate_ground_station_pairs(num_rings, num_sats_per_ring);
    int R_total = requirements.size();

    // --- GS SAT data ---
    vector<vector<GS_Sat_Data>> gs_sat_data(num_gs, vector<GS_Sat_Data>(num_sat));

    
    // --- Sat Loc ---
    
    double min_dis = numeric_limits<double>::max();

    for (int i = 0; i < num_sat; i++) {
        vector<double> sat_pos = satellites[i];
        // cout << "now on " << i << '\n';
        for (int j = 0; j < num_gs; j++) {
            vector<double> gs_pos = ground_stations[j];
            
            double tmp_dis = norm(subtract(sat_pos, gs_pos));
            double tmp_ang = compute_elevation_angle(gs_pos, sat_pos, true); // ���� (degree)
            double tmp_h_atm = ground_to_atmosphere_distance(gs_pos, sat_pos);
            double fid0 = 0.85;  // �T�w��
            // double height = norm(sat_pos) - R_earth_val;
            double tmp_fid = fidelity_Fij(tmp_ang, 20.0, eta_Tot(tmp_dis, tmp_h_atm), 1.0, fid0);
            // double tmp_fid = 1;
            double tmp_gen_rate = gen_rate(tmp_dis);
            if (tmp_dis < min_dis)
                min_dis = tmp_dis;
            gs_sat_data[j][i] = { static_cast<int>(tmp_dis), static_cast<int>(tmp_ang),
                                  static_cast<int>(tmp_h_atm), tmp_fid, tmp_gen_rate };
        }
    }

    cerr << "get all gs sat data\n";

    ofstream fout(output_filename);
    assert(fout);

    // fout << "S G R\n";
    fout << satellites.size() << " " << ground_stations.size() << " " << requirements.size() << "\n\n";

    // fout << "Satellite location\n";
    for (auto &sat : satellites) {
        fout << sat[0] << " " << sat[1] << " " << sat[2] << "\n";
    }
    // fout << "\nGS location\n";
    for (auto &gs : ground_stations) {
        fout << gs[0] << " " << gs[1] << " " << gs[2] << "\n";
    }
    // fout << "\n Request\n";
    for (auto &p : requirements) {
        fout << p.first << " " << p.second << "\n";
    }
    // fout << "\nData dis, ang, h_atm, fid, gen_rate\n";
    for (int i = 0; i < num_sat; i++) {
        for (int j = 0; j < num_gs; j++) {
            GS_Sat_Data data = gs_sat_data[j][i];
            fout << data.dis << " " << data.ang << " " << data.h_atm << " "
                 << data.fid << " " << data.gen_rate << "\n";
        }
    }
    fout.close();

    cout << "finish all output to dataset.txt" << endl;
    return 0;
}
