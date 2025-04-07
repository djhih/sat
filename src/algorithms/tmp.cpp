#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <utility>

using namespace std;

// 自訂 pair<int, int> 的 hash 結構
struct PairHash {
    size_t operator()(const pair<int, int>& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};

// 定義地心座標系 (ECEF)
struct ECEF {
    double x, y, z;
};

// 通道模型功能，計算 fidelity 與 transmittance
namespace QuantumChannel {
    inline double computeFidelity(double Fi0, double nj, double nij) {
        return 0.25 * (1 + (4 * Fi0 - 1) / pow(1 + nj / nij, 2));
    }

    inline double computeTransmittance(double dT, double dR, double lambda, double alpha, double h) {
        double pi = M_PI;
        double term1 = pi * pow(dT / 2, 2);
        double term2 = pi * pow(dR / 2, 2);
        double denominator = pow(lambda, 2);
        return (term1 * term2 / denominator) * exp(-alpha * h);
    }
}

// 衛星類別，封裝了與衛星相關的屬性與方法
class Satellite {
public:
    ECEF location;
    // 使用 map 來記錄衛星能服務的 ground station (GS) 與 GS pair，
    // key 可以是 GS 的 id 或 GS pair 的識別碼，value 為相關評估指標
    unordered_map<int, double> gsService;
    unordered_map<int, double> gsPairService;
    int servedReqId = -1;

    // 計算與某個 ground station 間的 fidelity (範例實作，可根據實驗模型完善)
    double countFidelity(const ECEF& gs) const {
        // TODO: 根據實際的距離與其他參數計算 fidelity
        return 0.0;
    }

    // 計算與某個 ground station 的仰角
    double countAngle(const ECEF& gs) const {
        // TODO: 根據位置計算衛星與 ground station 間的仰角
        return 0.0;
    }

    // 計算大氣衰減影響 (示意)
    double countAtmosphericLoss(const ECEF& gs) const {
        // TODO: 實作計算大氣損耗
        return 0.0;
    }

    // 根據設定的角度閾值檢查是否滿足條件
    bool checkAngle(double threshold) const {
        // TODO: 使用 countAngle() 並與 threshold 比較
        return true;
    }

    // 根據 fidelity 閾值檢查是否滿足條件
    bool checkFidelity(double threshold) const {
        // TODO: 使用 countFidelity() 並與 threshold 比較
        return true;
    }

    // 檢查同時滿足兩個 ground station 的條件
    bool checkBothGroundStations(double angleThreshold, double fidelityThreshold) const {
        return checkAngle(angleThreshold) && checkFidelity(fidelityThreshold);
    }
};

// GroundStation 類別
class GroundStation {
public:
    int id;
    ECEF location;
};

// Requirement 類別，代表一組 ground station pair 的需求
class Requirement {
public:
    int id;
    GroundStation gs1;
    GroundStation gs2;
    int servedBySatId = -1;
    double generationRate = 0.0; // 生成速率

    // 根據衛星提供的生成速率更新該需求的最小生成速率
    void updateGenerationRate(int satId,
        const unordered_map<pair<int, int>, double, PairHash>& rateGsSat,
        unordered_map<pair<int, int>, double, PairHash>& rateGsPairSat) {
        pair<int, int> key1 = {gs1.id, satId};
        pair<int, int> key2 = {gs2.id, satId};
        auto it1 = rateGsSat.find(key1);
        auto it2 = rateGsSat.find(key2);
        if (it1 != rateGsSat.end() && it2 != rateGsSat.end()) {
            generationRate = min(it1->second, it2->second);
            rateGsPairSat[{id, satId}] = generationRate;
        }
    }
};

// 模擬實驗的核心類別，管理整個衛星網路、ground station 以及需求
class SatelliteNetworkSimulation {
private:
    vector<Satellite> satellites;
    vector<GroundStation> groundStations;
    vector<Requirement> requirements;

    // 封裝衛星與 ground station 之間各種參數的資料結構
    unordered_map<pair<int, int>, double, PairHash> rateGsSat;
    unordered_map<pair<int, int>, double, PairHash> rateGsPairSat;
    unordered_map<pair<int, int>, double, PairHash> distanceGsSat;
    unordered_map<pair<int, int>, double, PairHash> angleGsSat;
    unordered_map<pair<int, int>, pair<double, double>, PairHash> fidelityGsSat;

    // 其他全域參數可以作為類別成員變數保存
    double fidelityThreshold = 0.5;
    double angleThreshold = 20.0;
    double alpha = 0.028125; // 大氣消光係數

public:
    // 構造函數，根據參數初始化各個容器
    SatelliteNetworkSimulation(int numSats, int numGS, int numReq) {
        satellites.resize(numSats);
        groundStations.resize(numGS);
        requirements.resize(numReq);
    }

    // 執行模擬實驗的主流程
    void runSimulation() {
        // TODO: 根據需求初始化衛星、GS 與需求參數
        // 例如：
        // - 計算衛星與 GS 之間的距離、仰角等參數，並存入 distanceGsSat 與 angleGsSat
        // - 根據 QuantumChannel::computeTransmittance() 計算傳輸效率
        // - 根據不同策略為每個需求分配最佳的衛星，更新 generationRate 等
        cout << "Running simulation with " << satellites.size() << " satellites, "
             << groundStations.size() << " ground stations, and "
             << requirements.size() << " requirements." << endl;
        
        // 示意：遍歷所有需求並更新生成速率 (需結合實際的 rateGsSat 資料)
        for (auto& req : requirements) {
            // 這裡選擇第一顆衛星作為範例
            int satId = 0;
            req.updateGenerationRate(satId, rateGsSat, rateGsPairSat);
        }
    }

    // 根據需要可加入其他輔助函數，例如初始化網路模型、策略選擇等
};

int main() {
    // 設定模擬參數
    int numSatellites = 100;
    int numGroundStations = 50;
    int numRequirements = 30;

    // 建立模擬實例並運行實驗
    SatelliteNetworkSimulation simulation(numSatellites, numGroundStations, numRequirements);
    simulation.runSimulation();

    return 0;
}
