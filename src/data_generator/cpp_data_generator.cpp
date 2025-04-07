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
// #include <direct.h>
#include <limits.h>
#include <string>
using namespace std;

const double PI = acos(-1.0);

// --- vector ---
double norm(const vector<double>& v) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double dot(const vector<double>& a, const vector<double>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

vector<double> subtract(const vector<double>& a, const vector<double>& b) {
    return { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
}

vector<double> scalar_divide(const vector<double>& v, double s) {
    return { v[0]/s, v[1]/s, v[2]/s };
}

// --- 基本物理與座標轉換 ---
// 依據球面座標 (R, theta, phi) 轉換成 Cartesian 座標 (x,y,z)
// 其中：theta 為環（固定極角），phi 為同一環中衛星分布角
vector<double> generate_point(double R, double theta, double phi) { 
    return { R * sin(phi) * cos(theta),
             R * sin(phi) * sin(theta),
             R * cos(phi) };
}

double get_R_earth() {
    return 6.371e6;  // 地球半徑 (m)
}

// --- 穿透率相關計算 ---
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

// --- 保真度與生成速率 ---
double fidelity_Fij(double theta_k1, double theta_e, double n_sg, double n_ij, double F0) {
    ofstream out("debug_fid.txt", std::ios::app);

    if ((theta_k1) >= theta_e) {
        out << "nsg " << n_sg << " n_ij " << n_ij << " F0 " << F0 << '\n';
        out << "Fid " << 0.25 * (1.0 + (4.0 * F0 - 1.0) / pow((1.0 + n_sg / n_ij), 2)) << '\n';
        return 0.25 * (1.0 + (4.0 * F0 - 1.0) / pow((1.0 + n_sg / n_ij), 2));
    } else {
        return 0.0;
    }
    out.close();
}

//! 
double gen_rate(double d) {
    // return 20.0;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(20, 50);
    return dis(gen);

    // double ratio = d / (6.371e6+500000);
    // return 50.0 * (ratio * ratio);
}

// --- 連結距離與最佳高度 ---
double link_distance(double d, double h) {
    double R = get_R_earth();
    return sqrt(4.0 * R * (R + h) * pow(sin(d / (4.0 * R)), 2) + h * h);
}

double opt_alt_arc(double d) { // not used
    auto obj = [d](double h) -> double {
        return -eta_Tot(link_distance(d, h), h);
    };

    double a = 0.0, b = 1500000.0;
    double tol = 1e-6;
    const int max_iter = 100;
    double gr = (sqrt(5.0) - 1) / 2.0;
    double c = b - gr * (b - a);
    double d_val = a + gr * (b - a);
    double fc = obj(c);
    double fd = obj(d_val);
    
    for (int i = 0; i < max_iter; i++) {
        if (fabs(b - a) < tol)
            break;
        if (fc < fd) {
            b = d_val;
            d_val = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = obj(c);
        } else {
            a = c;
            c = d_val;
            fc = fd;
            d_val = a + gr * (b - a);
            fd = obj(d_val);
        }
    }
    return (a + b) / 2.0;
}

double generate_link_constraint(double d) { // not used 2
    double h_opt = opt_alt_arc(d);
    return link_distance(d, h_opt);
}

// --- 仰角 ---
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

double ground_to_atmosphere_distance(const vector<double>& ground_ecef, const vector<double>& sat_ecef,
                                     double R = 6371000, double atmosphere_thickness = 10000) { 
    vector<double> p0 = ground_ecef;
    vector<double> v = subtract(sat_ecef, ground_ecef);
    double R_atm = R + atmosphere_thickness;
    double a = dot(v, v);
    double b = 2 * dot(p0, v);
    double c = dot(p0, p0) - R_atm * R_atm;
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        cerr << "無交點，請檢查輸入參數" << endl;
        return -1;
    }
    double sqrt_disc = sqrt(discriminant);
    double t1 = (-b + sqrt_disc) / (2 * a);
    double t2 = (-b - sqrt_disc) / (2 * a);
    double t_val = (t1 > 0) ? t1 : t2;
    double distance = t_val * norm(v);
    return distance;
}

// 依據環數與每環衛星數，利用球面座標產生衛星位置
// 與原 generate_network 中不同，此版本避免 k=0 產生重複點，改用 (k+1)*2*pi/num_sats 作為衛星在環內分布角
void generate_satellite_positions(int num_rings, int num_sats, double R, double h, double offset_angle,
                                  vector<vector<vector<double>>>& S,
                                  vector<vector<double>>& S_axis) {
    S.resize(num_rings);
    S_axis.resize(num_rings);
    for (int m = 0; m < num_rings; m++) {
        S[m].resize(num_sats);
        double theta_ring = m * PI / num_rings;
        for (int k = 0; k < num_sats; k++) {
            double phi = (k + 1) * 2 * PI / num_sats + offset_angle;
            S[m][k] = generate_point(R + h, theta_ring, phi);
        }
        S_axis[m] = { -sin(theta_ring), cos(theta_ring), 0.0 };
    }
}

void generate_network(int num_rings, int num_sats, double R, double h,
                      vector<vector<vector<double>>>& S,
                      vector<vector<vector<double>>>& G,
                      vector<vector<double>>& S_axis,
                      double offset_angle = 0.0) {
    // 產生衛星與旋轉軸（利用新函式）
    generate_satellite_positions(num_rings, num_sats, R, h, offset_angle, S, S_axis);
    
    G.resize(num_rings);
    for (int m = 0; m < num_rings; m++) {
        G[m].resize(num_sats);
        double theta_gs = m * PI / num_rings + PI / (2 * num_rings);
        for (int k = 0; k < num_sats; k++) {
            double phi = k * 2 * PI / num_sats + offset_angle;
            G[m][k] = generate_point(R, theta_gs, phi);
        }
    }
}

// not uesd
vector<pair<int, int>> generate_ground_station_pairs(int num_rings, int num_sats) {
    vector<pair<int, int>> pairs;
    for (int m = 0; m < num_rings; m++) {
        for (int k = 0; k < num_sats; k++) {
            int gs1 = m * num_sats + k;
            int gs2;
            if (m == num_rings - 1)
                gs2 = 0 * num_sats + ((k + 1) % num_sats);
            else
                gs2 = (m + 1) * num_sats + k;
            pairs.push_back({ gs1, gs2 });
        }
    }
    return pairs;
}

string gs_filename = "dataset/code/output/gs_loc.txt";
void input_GS_data(vector<vector<double>>& ret) {
    ifstream in(gs_filename);
    char buffer[1000];
    assert(in);  

    set<vector<double>> gs_sat;
    int GS;
    in >> GS;

    if (GS <= 0) { 
        cerr << "Error: Invalid number of ground stations!\n";
        assert(0);  
    }

    for (int i = 0; i < GS; i++) {
        vector<double> tmp(3);  
        in >> tmp[0] >> tmp[1] >> tmp[2];
        gs_sat.insert(tmp);
    }

    ret = vector<vector<double>>(gs_sat.begin(), gs_sat.end());

    cerr << "DEBUG: input_GS_data()\n"; 
}

vector<pair<int, int>> random_ground_satation_pair(int G_total, int GS) {
    vector<pair<int, int>> result;
    set<pair<int, int>> res;
    
    if (G_total < 2) {
        cerr << "Error: G_total must be at least 2!\n";
        return result; // 空的 vector
    }

    // 建立隨機數產生器
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, GS - 1);

    while(res.size() < G_total){
        int a = dis(gen);
        int b = dis(gen);
        while(a == b || res.count({b, a})){
            b = dis(gen);
        }
        res.insert({a, b});
    }

    result = vector<pair<int, int>>(res.begin(), res.end());

    return result;  // 確保有 return 一個 vector
}

vector<pair<int, int>> get_all_ground_station_pairs(int GS) {
    vector<pair<int, int>> pairs;
    for (int i = 0; i < GS; i++) {
        for (int j = i + 1; j < GS; j++) {
            pairs.push_back({i, j});
        }
    }
    return pairs;
}

struct GS_Sat_Data {
    int dis;
    int ang;
    int h_atm;
    double fid;
    double gen_rate;
};

int main(int argc, char* argv[]) {
    ofstream out("debug_fid.txt");
    out.clear();
    out.close();
    int num_rings = 20;            // 環數
    int num_sats_per_ring = 20;      // 每環衛星數
    double R_earth_val = get_R_earth();   // 地球半徑 (m)
    double h = 500e3;             // 衛星高度 (m)
    
    string output_filename = "dataset/raw/dataset.txt";

    if(argc == 5){
        num_rings = stoi(argv[1]);
        num_sats_per_ring = stoi(argv[2]);
        gs_filename = argv[3];
        output_filename = argv[4];
    } 


    vector<vector<vector<double>>> S_dict;
    vector<vector<vector<double>>> G_dict;
    vector<vector<double>> S_axis;
    generate_network(num_rings, num_sats_per_ring, R_earth_val, h, S_dict, G_dict, S_axis);

    vector<vector<double>> satellites;
    for (int m = 0; m < num_rings; m++) {
        for (int k = 0; k < num_sats_per_ring; k++) {
            satellites.push_back(S_dict[m][k]);
        }
    }
    
    // 處理重複點
    set<vector<double>> S_set(satellites.begin(), satellites.end());
    satellites = vector<vector<double>>(S_set.begin(), S_set.end());
    vector<vector<double>> sat_parameters = satellites;
    int S_total = satellites.size();
    assert(satellites.size() == sat_parameters.size());
    assert(S_total == sat_parameters.size());
    
    // --- GS loc ---
    vector<vector<double>> ground_stations;
    input_GS_data(ground_stations);
    int G_total = ground_stations.size();
    assert(G_total != 0);
    assert(G_total == ground_stations.size());

    // --- GSP ---
    // vector<pair<int, int>> requirements = random_ground_satation_pair(10000, G_total);
    vector<pair<int, int>> requirements = get_all_ground_station_pairs(G_total);

    // vector<pair<int, int>> requirements = generate_ground_station_pairs(num_rings, num_sats_per_ring);
    int R_total = requirements.size();


    // --- GS SAT data ---
    vector<vector<GS_Sat_Data>> gs_sat_data(G_total, vector<GS_Sat_Data>(S_total));

    
    // --- Sat Loc ---
    
    double min_dis = numeric_limits<double>::max();

    for (int i = 0; i < S_total; i++) {
        vector<double> sat_pos = satellites[i];
        for (int j = 0; j < G_total; j++) {
            vector<double> gs_pos = ground_stations[j];
            double tmp_dis = norm(subtract(sat_pos, gs_pos));
            double tmp_ang = compute_elevation_angle(gs_pos, sat_pos, true); // 仰角 (degree)
            double tmp_h_atm = ground_to_atmosphere_distance(gs_pos, sat_pos, R_earth_val, 10000);
            double fid0 = 0.85;  // 固定值
            double height = norm(sat_pos) - R_earth_val;
            double tmp_fid = fidelity_Fij(tmp_ang, 20.0, eta_Tot(tmp_dis, height), 1.0, fid0);
            double tmp_gen_rate = gen_rate(tmp_dis);
            if (tmp_dis < min_dis)
                min_dis = tmp_dis;
            gs_sat_data[j][i] = { static_cast<int>(tmp_dis), static_cast<int>(tmp_ang),
                                  static_cast<int>(tmp_h_atm), tmp_fid, tmp_gen_rate };
        }
    }

    
    

    ofstream fout(output_filename);
    assert(fout);

    // fout << "S G R\n";
    fout << sat_parameters.size() << " " << ground_stations.size() << " " << requirements.size() << "\n\n";

    // fout << "Satellite location\n";
    for (auto &sat : sat_parameters) {
        fout << sat[0] << " " << sat[1] << " " << sat[2] << "\n";
    }
    // // fout << "\nGS location\n";
    for (auto &gs : ground_stations) {
        fout << gs[0] << " " << gs[1] << " " << gs[2] << "\n";
    }
    // fout << "\n Request\n";
    for (auto &p : requirements) {
        fout << p.first << " " << p.second << "\n";
    }
    // fout << "\nData dis, ang, h_atm, fid, gen_rate\n";
    for (int i = 0; i < S_total; i++) {
        for (int j = 0; j < G_total; j++) {
            GS_Sat_Data data = gs_sat_data[j][i];
            fout << data.dis << " " << data.ang << " " << data.h_atm << " "
                 << data.fid << " " << data.gen_rate << "\n";
            // if(data.fid != 0){
            //     cout << data.dis << " " << data.ang << " " << data.h_atm << " " << data.fid << " " << data.gen_rate << "\n";
            // }
        }
    }
    fout.close();

    // cout << "output to dataset.txt" << endl;
    return 0;
}
