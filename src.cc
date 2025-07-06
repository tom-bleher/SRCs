#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <omp.h>
#include <memory>

using namespace std;

// Global parameters
double pairProb = 0.2;
double Rsrc = 0.5;

// Pair accumulation arrays (retain original functionality)
double np_array[4][8][4][8] = {};
double pp_array[4][8][4][8] = {};
double nn_array[4][8][4][8] = {};

// Optimized data structures
struct alignas(32) NucleonPosition {
    float x, y, z;
    float probability;
    float transparency;
    uint16_t type;           // 0 = proton, 1 = neutron
    uint16_t shell_index;    // index into shell states
    uint16_t n_r, l;         // quantum numbers for fast access
};

// Global transparency function to match original
double getTransparency(int nucl, double x, double y, double z, double rmax, double sigmaNN) {
    // (x,y,z) = initial cell position
    double transp= 0.0;
    double dStep = 0.01;
    double rPrime;
    do { // integral rho sigma dz
        rPrime=std::sqrt(x*x+y*y+z*z);
        double dens = 0.0;
        if(nucl == 1){
            dens = 0.0922/(1+exp((rPrime-2.861)/0.52));
        }
        if(nucl == 2){
            dens = 0.1027/(1+exp((rPrime-3.75)/0.52));
        }
        if(nucl == 3){
            dens = 0.1095/(1+exp((rPrime-4.781)/0.52));
        }
        if(nucl == 4){
            dens = 0.1156/(1+exp((rPrime-7.423)/0.52));
        }
        if(nucl == 5){
            dens = 0.1066/(1+exp((rPrime-4.274)/0.52));
        }
        if(nucl == 6){
            dens = 0.1081/(1+exp((rPrime-4.543)/0.52));
        }
        if(nucl == 7){
            dens = 0.1090/(1+exp((rPrime-4.724)/0.52));
        }
        transp = transp + (sigmaNN*dens*dStep);
        z = z + dStep;
    } while (rPrime<rmax);
    transp = exp(-transp);
    return transp;
}

// Original wavefunction calculation functions (exact match to original)
double assoc_laguerre(int n, double alpha, double x) {
    double sum = 0.0;
    for (int k = 0; k <= n; ++k) {
        double num = std::tgamma(n + alpha + 1);
        double den = std::tgamma(n - k + 1) * std::tgamma(alpha + k + 1);
        double binom = num / den;
        sum += std::pow(-x, k) * binom / std::tgamma(k + 1);
    }
    return sum;
}

double normalization(int n_r, int l, double b) {
    double num = 2.0 * std::tgamma(n_r + 1);
    double den = std::pow(b, 3) * std::tgamma(n_r + l + 1.5);
    return std::sqrt(num / den);
}

double R_radial(int n_r, int l, double b, double r) {
    double rho = r / b;
    double rho2 = rho * rho;
    double alpha = l + 0.5;
    double L = assoc_laguerre(n_r, alpha, rho2);
    double N = normalization(n_r, l, b);
    return N * std::pow(rho, l) * L * std::exp(-rho2 / 2.0);
}

// Structure for a shell state
struct ShellState {
    int n_r;
    int l;
};

// Cross section calculation
double getReep(double np, double pp, double nn, double Tnp_p, double Tpp1, double Tpp2, double Tnp_n, double Tnn1, double Tnn2) {
    const double P = 0.05;
    const double sigmap = 2.4;
    const double sigman = 1.0;
    
    return (Tnp_p*np+(Tpp1+Tpp2)*pp)*sigmap + (Tnp_n*np+(Tnn1+Tnn2)*nn)*sigman*P;
}

int main(int argc, char* argv[]) {
    if(argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <nucleus_type> <Rsrc> <pairProb>\n";
        return 1;
    }
    
    pairProb = std::stod(argv[3]);
    int nuclType = std::stoi(argv[1]);
    Rsrc = std::stod(argv[2]);
    const double Rsrc_squared = Rsrc * Rsrc;
    
    // Optimized parameters
    const double r_max = 10.0;
    const double spacing = 0.3;
    const double sigmaNN = 3.5;
    const double LimitAngular = 100.0;
    
    double A = 208.0;
    
    // Nuclear mass setup
    const double mass_numbers[] = {2, 12, 27, 56, 208, 40, 48, 54, 9, 63, 197};
    if (nuclType >= 0 && nuclType <= 10) {
        A = mass_numbers[nuclType];
    }
    
    // Oscillator parameters
    const double hbarc = 197.3269804;
    const double m_n = 939.0;
    double hbar_omega = 41.0 * std::pow(A, -1.0/3.0);
    double b = hbarc / std::sqrt(m_n * hbar_omega);

    // Initialize nuclear configurations
    std::vector<ShellState> protons, neutrons;
    double b_p, b_n;
    
    // Complete nuclear configuration setup
    if(nuclType == 0) { // 2H
        protons = {{0,0}};
        neutrons = {{0,0}};
        
        double hbar_omega_p = 41.0 * std::pow(1, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(1, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    else if(nuclType == 1) { // 12C
        protons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}, {0,1}};
        neutrons = protons;
        
        double hbar_omega_p = 41.0 * std::pow(6, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(6, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    else if(nuclType == 2) { // 27Al
        protons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2}};
        neutrons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2}};
        
        double hbar_omega_p = 41.0 * std::pow(13, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(14, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    else if(nuclType == 3) { // 56Fe
        protons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{1,0},{1,0},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3}};
        neutrons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{1,0},{1,0},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{1,1},{1,1}};
        
        double hbar_omega_p = 41.0 * std::pow(26, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(30, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    else if(nuclType == 4) { // 208Pb
        protons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{1,0},{1,0},{0,2},{0,2},{0,2},{0,2},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{1,1},{1,1},{1,1},{1,1},{0,3},{0,3},{1,1},{1,1},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{1,2},{1,2},{1,2},{1,2},{1,2},{1,2},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{2,0},{2,0},{1,2},{1,2},{1,2},{1,2},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5}};
        neutrons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{1,0},{1,0},{0,2},{0,2},{0,2},{0,2},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{1,1},{1,1},{1,1},{1,1},{0,3},{0,3},{1,1},{1,1},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{1,2},{1,2},{1,2},{1,2},{1,2},{1,2},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{2,0},{2,0},{1,2},{1,2},{1,2},{1,2},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{1,3},{1,3},{1,3},{1,3},{1,3},{1,3},{1,3},{1,3},{2,1},{2,1},{2,1},{2,1},{1,3},{1,3},{1,3},{1,3},{1,3},{1,3},{2,1},{2,1},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{2,4},{2,4},{2,4},{2,4},{2,4},{2,4},{2,4},{2,4},{2,4},{2,4}};
        
        double hbar_omega_p = 41.0 * std::pow(82, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(126, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    else if(nuclType == 5) { // 40Ca
        protons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{1,0},{1,0}};
        neutrons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{1,0},{1,0}};
        
        double hbar_omega_p = 41.0 * std::pow(20, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(20, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    else if(nuclType == 6) { // 48Ca
        protons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{1,0},{1,0}};
        neutrons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{1,0},{1,0},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3}};
        
        double hbar_omega_p = 41.0 * std::pow(20, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(28, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    else if(nuclType == 7) { // 54Fe
        protons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{1,0},{1,0},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3}};
        neutrons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{1,0},{1,0},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3}};
        
        double hbar_omega_p = 41.0 * std::pow(26, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(28, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    else if(nuclType == 8) { // 9Be
        protons = {{0,0}, {0,0}, {0,1}, {0,1}};
        neutrons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}};
        
        double hbar_omega_p = 41.0 * std::pow(4, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(5, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    else if(nuclType == 9) { // 63Cu
        protons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {1,0}, {1,0}, {0,2}, {0,2}, {0,2}, {0,2}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {1,1}};
        neutrons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {1,0}, {1,0}, {0,2}, {0,2}, {0,2}, {0,2}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {1,1}, {1,1}, {1,1}, {1,1}, {0,3}, {0,3}};
        
        double hbar_omega_p = 41.0 * std::pow(29, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(34, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    else if(nuclType == 10) { // 197Au
        protons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {1,0}, {1,0}, {0,2}, {0,2}, {0,2}, {0,2}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {1,1}, {1,1}, {1,1}, {1,1}, {0,3}, {0,3}, {1,1}, {1,1}, {2,0}, {2,0}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {1,3}, {1,3}, {1,3}, {1,3}, {1,3}, {1,3}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,2}, {2,2}, {2,2}};
        neutrons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {1,0}, {1,0}, {0,2}, {0,2}, {0,2}, {0,2}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {0,3}, {1,1}, {1,1}, {1,1}, {1,1}, {0,3}, {0,3}, {1,1}, {1,1}, {2,0}, {2,0}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {1,3}, {1,3}, {1,3}, {1,3}, {1,3}, {1,3}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,2}, {2,2}, {2,2}, {2,2}, {2,2}, {2,2}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {2,3}, {2,3}, {2,3}, {2,3}, {3,0}, {3,0}, {2,2}, {2,2}, {2,2}, {2,2}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {2,3}, {2,3}, {2,3}, {2,3}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}};
        
        double hbar_omega_p = 41.0 * std::pow(79, -1.0/3.0);
        double hbar_omega_n = 41.0 * std::pow(118, -1.0/3.0);
        b_p = hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n = hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    
    // Generate spatial grid more efficiently
    std::vector<std::array<double, 3>> centers;
    int grid_size = static_cast<int>(2 * r_max / spacing) + 1;
    centers.reserve(grid_size * grid_size * grid_size);
    
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            for (int k = 0; k < grid_size; ++k) {
                double x = -r_max + i * spacing;
                double y = -r_max + j * spacing;
                double z = -r_max + k * spacing;
                centers.push_back({x, y, z});
            }
        }
    }

    const double cellVol = spacing * spacing * spacing;
    const double accuracy = 1e-3;
    const double total_volume = (2.0 * r_max) * (2.0 * r_max) * (2.0 * r_max);
    const double N_nucleons = A * total_volume / (spacing * spacing * spacing);
    const double N_pairs = N_nucleons * (N_nucleons - 1) / 2;
    const double prob_threshold = std::sqrt(accuracy / N_pairs);
    
    std::cout << "Probability threshold: " << prob_threshold << std::endl;
    
    // Pre-compute transparency for every center (expensive integral done once)
    std::vector<double> centerTransparency;
    centerTransparency.resize(centers.size());
    
    #pragma omp parallel for schedule(static)
    for (size_t c_idx = 0; c_idx < centers.size(); ++c_idx) {
        const auto &c = centers[c_idx];
        centerTransparency[c_idx] = getTransparency(nuclType, c[0], c[1], c[2], r_max, sigmaNN);
    }

    // Create optimized nucleon data structure
    std::vector<NucleonPosition> nucleons;
    nucleons.reserve(protons.size() * centers.size() / 10);  // Estimate after filtering
    
    std::cout << "Creating nucleon probability distribution..." << std::endl;
    
    // Parallelize nucleon creation with better memory management
    #pragma omp parallel
    {
        std::vector<NucleonPosition> local_nucleons;
        local_nucleons.reserve(1000);
        
        #pragma omp for schedule(dynamic, 1)
        for (size_t p_idx = 0; p_idx < protons.size(); ++p_idx) {
            const auto& p = protons[p_idx];
            
            for (size_t c_idx = 0; c_idx < centers.size(); ++c_idx) {
                const auto& c = centers[c_idx];
                double r = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
                
                if (r > r_max) continue;
                
                double Rp = R_radial(p.n_r, p.l, b_p, r);
                double wf_sq = Rp * Rp;
                double prob = wf_sq / (4.0 * M_PI) * cellVol;
                
                if (prob > prob_threshold) {
                    double transp = centerTransparency[c_idx];
                    
                    NucleonPosition nucleon;
                    nucleon.x = static_cast<float>(c[0]);
                    nucleon.y = static_cast<float>(c[1]);
                    nucleon.z = static_cast<float>(c[2]);
                    nucleon.probability = static_cast<float>(prob);
                    nucleon.transparency = static_cast<float>(transp);
                    nucleon.type = 0;  // proton
                    nucleon.shell_index = static_cast<uint16_t>(p_idx);
                    nucleon.n_r = static_cast<uint16_t>(p.n_r);
                    nucleon.l = static_cast<uint16_t>(p.l);
                    
                    local_nucleons.push_back(nucleon);
                }
            }
        }
        
        #pragma omp for schedule(dynamic, 1)
        for (size_t n_idx = 0; n_idx < neutrons.size(); ++n_idx) {
            const auto& n = neutrons[n_idx];
            
            for (size_t c_idx = 0; c_idx < centers.size(); ++c_idx) {
                const auto& c = centers[c_idx];
                double r = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
                
                if (r > r_max) continue;
                
                double Rn = R_radial(n.n_r, n.l, b_n, r);
                double wf_sq = Rn * Rn;
                double prob = wf_sq / (4.0 * M_PI) * cellVol;
                
                if (prob > prob_threshold) {
                    double transp = centerTransparency[c_idx];
                    
                    NucleonPosition nucleon;
                    nucleon.x = static_cast<float>(c[0]);
                    nucleon.y = static_cast<float>(c[1]);
                    nucleon.z = static_cast<float>(c[2]);
                    nucleon.probability = static_cast<float>(prob);
                    nucleon.transparency = static_cast<float>(transp);
                    nucleon.type = 1;  // neutron
                    nucleon.shell_index = static_cast<uint16_t>(n_idx);
                    nucleon.n_r = static_cast<uint16_t>(n.n_r);
                    nucleon.l = static_cast<uint16_t>(n.l);
                    
                    local_nucleons.push_back(nucleon);
                }
            }
        }
        
        // Combine thread results
        #pragma omp critical
        {
            nucleons.insert(nucleons.end(), local_nucleons.begin(), local_nucleons.end());
        }
    }
    
    std::cout << "Created " << nucleons.size() << " nucleon-position combinations" << std::endl;

    // Start pair processing
    std::cout << "Computing pair probabilities..." << std::endl;
    
    // Results storage
    std::vector<double> thread_P_np(omp_get_max_threads(), 0.0);
    std::vector<double> thread_P_pp(omp_get_max_threads(), 0.0);
    std::vector<double> thread_P_nn(omp_get_max_threads(), 0.0);
    std::vector<double> thread_Reep(omp_get_max_threads(), 0.0);
    
    // Simplified pair arrays (reduce memory usage)
    std::array<std::array<std::array<std::array<double, 8>, 4>, 8>, 4> np_pairs{};
    std::array<std::array<std::array<std::array<double, 8>, 4>, 8>, 4> pp_pairs{};
    std::array<std::array<std::array<std::array<double, 8>, 4>, 8>, 4> nn_pairs{};
    
    // Use original data structures for exact compatibility
    std::vector<std::array<double,3>> nucleon_centers;
    nucleon_centers.reserve(nucleons.size());
    for (const auto& n : nucleons) {
        nucleon_centers.push_back({static_cast<double>(n.x), static_cast<double>(n.y), static_cast<double>(n.z)});
    }
    
    // Exact same nested loop structure as original
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        size_t pairs_processed = 0;
        
        #pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < nucleons.size(); ++i) {
            const auto& n1 = nucleons[i];
            const auto& c1 = nucleon_centers[i];
            
            for (size_t j = i + 1; j < nucleons.size(); ++j) {
                const auto& n2 = nucleons[j];
                const auto& c2 = nucleon_centers[j];
                
                // Calculate distance between centers (exact same as original)
                double dx = c1[0] - c2[0];
                double dy = c1[1] - c2[1];
                double dz = c1[2] - c2[2];
                double dist_squared = dx*dx + dy*dy + dz*dz;
                
                if (dist_squared >= Rsrc_squared) continue;
                
                // Skip same nucleon
                if (n1.type == n2.type && n1.shell_index == n2.shell_index) continue;
                
                pairs_processed++;
                if (pairs_processed % 100000 == 0 && thread_id == 0) {
                    std::cout << "Processed " << pairs_processed << " pairs (thread " << thread_id << ")\r" << std::flush;
                }
                
                // Calculate pair probability
                double pair_prob = static_cast<double>(n1.probability * n2.probability);
                
                // Process pair based on type
                if (n1.type != n2.type) {
                    // n-p pair
                    if (std::abs(static_cast<int>(n1.l) - static_cast<int>(n2.l)) <= LimitAngular) {
                        thread_P_np[thread_id] += pair_prob;
                        
                        double Tnp_p = (n1.type == 0) ? n1.transparency : n2.transparency;
                        double Tnp_n = (n1.type == 1) ? n1.transparency : n2.transparency;
                        
                        thread_Reep[thread_id] += getReep(pair_prob, 0, 0, Tnp_p, 0, 0, Tnp_n, 0, 0);

                        // Accumulate into global np array
                        #pragma omp atomic
                        np_array[(int)n2.n_r][(int)n2.l][(int)n1.n_r][(int)n1.l] += pair_prob;
                    }
                } else {
                    // Same type pair
                    if (std::abs(static_cast<int>(n1.l) - static_cast<int>(n2.l)) <= LimitAngular) {
                        double corr_prob = pair_prob * pairProb;
                        
                        if (n1.type == 0) {
                            thread_P_pp[thread_id] += corr_prob;
                            thread_Reep[thread_id] += getReep(0, corr_prob, 0, 0, n1.transparency, n2.transparency, 0, 0, 0);

                            // Accumulate into global pp array
                            #pragma omp atomic
                            pp_array[(int)n1.n_r][(int)n1.l][(int)n2.n_r][(int)n2.l] += corr_prob;
                        } else {
                            thread_P_nn[thread_id] += corr_prob;
                            thread_Reep[thread_id] += getReep(0, 0, corr_prob, 0, 0, 0, 0, n1.transparency, n2.transparency);

                            // Accumulate into global nn array
                            #pragma omp atomic
                            nn_array[(int)n1.n_r][(int)n1.l][(int)n2.n_r][(int)n2.l] += corr_prob;
                        }
                    }
                }
            }
        }
    }
    
    // Combine results
    double P_np = 0.0, P_pp = 0.0, P_nn = 0.0, Reep = 0.0;
    for (int t = 0; t < omp_get_max_threads(); ++t) {
        P_np += thread_P_np[t];
        P_pp += thread_P_pp[t];
        P_nn += thread_P_nn[t];
        Reep += thread_Reep[t];
    }
    
    // Calculate total probabilities
    double Pp = 0.0, Pn = 0.0;
    for (const auto& n : nucleons) {
        if (n.type == 0) Pp += n.probability;
        else Pn += n.probability;
    }

    // Output results
    std::cout << "\nTotal close-proximity probabilities (small spheres):\n";
    std::cout << "  #np: " << P_np << std::endl;
    std::cout << "  #pp: " << P_pp << std::endl;
    std::cout << "  #nn: " << P_nn << std::endl;
    std::cout << "  R(e,e'p): " << Reep << std::endl;
    
    std::cout << "rmax: " << r_max << "\n";
    std::cout << " SRC fraction = " << (P_np+P_pp+P_nn)/A << std::endl;
    std::cout << " pp/np = " << P_pp/P_np << "\t nn/np = " << P_nn/P_np << std::endl;
    std::cout << " Pp = " << Pp << " Pn " << Pn << std::endl;
    std::cout << "spacing: " << spacing << "\n";
    std::cout << " Rsrc = " << Rsrc << std::endl;
    std::cout << " pairProb = " << pairProb << std::endl;
    
    return 0;
}

