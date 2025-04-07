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
#include <string>
#include <stdexcept>
using namespace std;

const double PI = acos(-1.0);

int num_sat = 20;
int num_gs = 30;
int num_req = 900;

double dot(const vector<double>& a, const vector<double>& b) {
    double dp = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        dp += a[i] * b[i];
    }
    return dp;
}

vector<double> subtract(const vector<double>& a, const vector<double>& b) {
    return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
}

vector<double> scalar_divide(const vector<double>& v, double s) {
    return { v[0] / s, v[1] / s, v[2] / s };
}

double norm(const vector<double>& v) {
    return sqrt(dot(v, v));
}

vector<double> normalize(const vector<double>& v) {
    double n = norm(v);
    if (n == 0) 
        throw runtime_error("Zero vector cannot be normalized");
    vector<double> u(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        u[i] = v[i] / n;
    }
    return u;
}

// (R, theta, phi) -> (x, y, z)
vector<double> generate_point(double R, double theta, double phi) { 
    return { R * sin(phi) * cos(theta),
             R * sin(phi) * sin(theta),
             R * cos(phi) };
}

double get_R_earth() {
    return 6.371e6; // 地球半徑 (m)
}

double eta_sg(double L, double lam = 810e-9, double w0 = 0.025, double rg = 0.75) {
    double LR = PI * w0 * w0 / lam;
    double wg = w0 * sqrt(1 + (L * L) / (LR * LR));
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

#ifdef DEBUG
// 若開啟 DEBUG，則用此輔助函式寫除錯資訊，避免重複開關檔案
void write_debug(const string& msg) {
    static ofstream debug("debug_fid.txt", ios::app);
    debug << msg;
}
#endif

double fidelity_Fij(double theta_k1, double theta_e, double n_sg, double n_ij, double F0) {
    if (theta_k1 >= theta_e) {
#ifdef DEBUG
        write_debug("nsg " + to_string(n_sg) + " n_ij " + to_string(n_ij) + " F0 " + to_string(F0) + "\n");
        double fidVal = 0.25 * (1.0 + (4.0 * F0 - 1.0) / pow((1.0 + n_sg / n_ij), 2));
        write_debug("Fid " + to_string(fidVal) + "\n");
        return fidVal;
#else
        return 0.25 * (1.0 + (4.0 * F0 - 1.0) / pow((1.0 + n_sg / n_ij), 2));
#endif
    } 
    return 0.0;
}

double gen_rate(double /*d*/) {
    static mt19937 gen(random_device{}());
    static uniform_int_distribution<> dis(20, 50);
    return dis(gen);
}

double link_distance(double d, double h) {
    double R = get_R_earth();
    return sqrt(4.0 * R * (R + h) * pow(sin(d / (4.0 * R)), 2) + h * h);
}

double compute_elevation_angle(const vector<double>& ground_ecef, const vector<double>& sat_ecef, bool degrees = false) {
    vector<double> diff = subtract(sat_ecef, ground_ecef);
    double norm_diff = norm(diff);
    vector<double> up = scalar_divide(ground_ecef, norm(ground_ecef));
    double dot_up = dot(diff, up);
    double elev_rad = asin(dot_up / norm_diff);
    return degrees ? elev_rad * 180.0 / PI : elev_rad;
}

struct GS_Sat_Data {
    int dis;
    int ang;
    int h_atm;
    double fid;
    double gen_rate;
};

// 全域容器
vector<vector<double>> satellites;
vector<vector<double>> ground_stations;

double ground_to_atmosphere_distance(const vector<double>& gs_ecef, const vector<double>& sat_ecef) {
    const double R = 6371000;
    const double atmosphere_thickness = 10000;
    vector<double> d = subtract(sat_ecef, gs_ecef);
    d = normalize(d);
    
    double R_atm = R + atmosphere_thickness;
    double gs_dot_d = dot(gs_ecef, d);
    double gs_norm_sq = dot(gs_ecef, gs_ecef);
    double A = 1.0, B = 2.0 * gs_dot_d, C = gs_norm_sq - R_atm * R_atm;
    double discriminant = B * B - 4 * A * C;
    if (discriminant < 0) {
        cerr << "gs ecef " << gs_ecef[0] << " " << gs_ecef[1] << " " << gs_ecef[2] << "\n";
        cerr << "sat ecef " << sat_ecef[0] << " " << sat_ecef[1] << " " << sat_ecef[2] << "\n";
        cerr << "gs norm " << gs_norm_sq << "\n";
        throw runtime_error("No intersection: 地面站必在大氣層內");
    }
    double sqrt_disc = sqrt(discriminant);
    double t1 = (-B - sqrt_disc) / (2 * A);
    double t2 = (-B + sqrt_disc) / (2 * A);
    double t = (t1 >= 0) ? t1 : t2;
    if (t < 0)
        throw runtime_error("射線與大氣層無正向交點");
    return t;
}


string sat_filename = "dataset/code/output/satellite_coordinates_selected.txt";
void input_sat_data(vector<vector<double>>& sats) {
    ifstream in(sat_filename);
    if (!in) throw runtime_error("無法開啟 " + sat_filename);
    in >> num_sat;
    sats.resize(num_sat, vector<double>(3));
    for (int i = 0; i < num_sat; i++) {
        in >> sats[i][0] >> sats[i][1] >> sats[i][2];
    }
}

string gs_filename = "dataset/code/output/gs_loc.txt";
void input_GS_data(vector<vector<double>>& gs_data) {
    ifstream in(gs_filename);
    if (!in) throw runtime_error("無法開啟 " + gs_filename);
    set<vector<double>> gs_set;
    in >> num_gs;
    for (int i = 0; i < num_gs; i++) {
        vector<double> tmp(3);
        in >> tmp[0] >> tmp[1] >> tmp[2];
        gs_set.insert(tmp);
    }
    gs_data.assign(gs_set.begin(), gs_set.end());
}


vector<pair<int, int>> random_ground_satation_pair(int num_pairs, int G_total) {
    vector<pair<int, int>> result;
    set<pair<int, int>> res;
    
    if (G_total < 2) {
        cerr << "Error: G_total must be at least 2!\n";
        return result;
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, num_gs - 1);

    set<pair<int, int>> pair_tried;
    while (res.size() < static_cast<size_t>(num_pairs) && pair_tried.size() < static_cast<size_t>(G_total * G_total)) {
        int a = dis(gen), b = dis(gen);
        while (a == b || res.count({b, a})) {
            b = dis(gen);
        }
        pair_tried.insert({a, b});
        // 檢查是否存在一顆衛星能同時讓兩個地面站的 fidelity > 0.5
        bool can_serve = false;
        for (const auto& sat : satellites) {
            double tmp_dis_a = norm(subtract(sat, ground_stations[a]));
            double tmp_dis_b = norm(subtract(sat, ground_stations[b]));
            double tmp_ang_a = compute_elevation_angle(ground_stations[a], sat, true);
            double tmp_ang_b = compute_elevation_angle(ground_stations[b], sat, true);
            double tmp_h_atm_a = ground_to_atmosphere_distance(ground_stations[a], sat);
            double tmp_h_atm_b = ground_to_atmosphere_distance(ground_stations[b], sat);
            double fid0 = 0.85;
            double tmp_fid_a = fidelity_Fij(tmp_ang_a, 20.0, eta_Tot(tmp_dis_a, tmp_h_atm_a), 1.0, fid0);
            double tmp_fid_b = fidelity_Fij(tmp_ang_b, 20.0, eta_Tot(tmp_dis_b, tmp_h_atm_b), 1.0, fid0);
            if (tmp_fid_a > 0.5 && tmp_fid_b > 0.5) {
                can_serve = true;
                break;
            }
        }
        res.insert({a, b});
    }
    result.assign(res.begin(), res.end());
    cerr << "試了 " << pair_tried.size() << " 組配對，得到 " << result.size() << " 組\n";
    return result;
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

int main(int argc, char* argv[]) {
    
    double R_earth_val = get_R_earth();   // 地球半徑 (m)
    double h = 500e3;                     // 衛星高度 (m)
    string output_filename = "dataset/raw/dataset.txt";

    if (argc == 5) {
        num_req = stoi(argv[1]);
        sat_filename = argv[2];
        gs_filename = argv[3];
        output_filename = argv[4];
    } 

    // --- 輸入衛星資料 ---
    input_sat_data(satellites);
    cerr << "取得 " << satellites.size() << " 筆衛星資料\n";

    // --- 輸入地面站資料 ---
    input_GS_data(ground_stations);
    cerr << "取得地面站資料\n";

    // --- 產生地面站配對 ---
    vector<pair<int, int>> requirements = random_ground_satation_pair(num_req, num_gs);
    cerr << "取得 " << requirements.size() << " 組地面站配對\n";

    // --- 建立 GS 與衛星間資料的二維矩陣 ---
    vector<vector<GS_Sat_Data>> gs_sat_data(num_gs, vector<GS_Sat_Data>(num_sat));

    double min_dis = numeric_limits<double>::max();

    // --- 計算各衛星與各地面站之數值 ---
    for (int i = 0; i < num_sat; i++) {
        const auto& sat_pos = satellites[i];
        for (int j = 0; j < num_gs; j++) {
            const auto& gs_pos = ground_stations[j];
            double tmp_dis = norm(subtract(sat_pos, gs_pos));
            double tmp_ang = compute_elevation_angle(gs_pos, sat_pos, true); // 仰角 (degree)
            double tmp_h_atm = ground_to_atmosphere_distance(gs_pos, sat_pos);
            double fid0 = 0.85;  // 固定基準值
            double tmp_fid = fidelity_Fij(tmp_ang, 20.0, eta_Tot(tmp_dis, tmp_h_atm), 1.0, fid0);
            double tmp_gen_rate = gen_rate(tmp_dis);
            if (tmp_dis < min_dis)
                min_dis = tmp_dis;
            gs_sat_data[j][i] = { static_cast<int>(tmp_dis), static_cast<int>(tmp_ang),
                                  static_cast<int>(tmp_h_atm), tmp_fid, tmp_gen_rate };
        }
    }
    cerr << "完成 GS 與衛星資料計算\n";

    // --- 將結果輸出至檔案 ---
    ofstream fout(output_filename);
    if (!fout) throw runtime_error("無法開啟輸出檔案 " + output_filename);

    fout << satellites.size() << " " << ground_stations.size() << " " << requirements.size() << "\n\n";
    for (const auto& sat : satellites)
        fout << sat[0] << " " << sat[1] << " " << sat[2] << "\n";
    for (const auto& gs : ground_stations)
        fout << gs[0] << " " << gs[1] << " " << gs[2] << "\n";
    for (const auto& p : requirements)
        fout << p.first << " " << p.second << "\n";
    for (int i = 0; i < num_sat; i++) {
        for (int j = 0; j < num_gs; j++) {
            const GS_Sat_Data& data = gs_sat_data[j][i];
            fout << data.dis << " " << data.ang << " " << data.h_atm << " "
                 << data.fid << " " << data.gen_rate << "\n";
        }
    }
    fout.close();

    cout << "完成所有資料輸出至 " << output_filename << endl;
    return 0;
}
