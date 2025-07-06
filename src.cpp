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

using namespace std;

double pairProb = 0.2;
double Rsrc = 0.5;

double np[4][8][4][8]; //n: n,l, p: n,l
double pp[4][8][4][8]; //
double nn[4][8][4][8];

struct TupleHash {
    size_t operator()(const std::tuple<int, int, int>& t) const {
        auto [x, y, z] = t;
        return std::hash<int>()(x) ^ (std::hash<int>()(y) << 1) ^ (std::hash<int>()(z) << 2);
    }
};

// Structure for a shell state (n_r, l)
struct ShellState {
    int n_r;
    int l;
};

// Structure to represent a nucleon at a specific position
struct Nucleon {
    int type; // 0 = proton, 1 = neutron
    int shell_index; // index into protons or neutrons vector
    size_t center_index; // which spatial center it's at
    double probability; // probability of being at this center
    double transparency; // transparency of nucleon from current position
};

double GenDensity(int nucl, double r){ // 1- carbon, 2- al, 3 - Fe56, 4 - Pb, 5 - Ca40, 6- Ca48 ,7 - Fe54
    double dens = 0.0;
    if(nucl == 1){
        dens = 0.0922/(1+exp((r-2.861)/0.52));
    }
    if(nucl == 2){
        dens = 0.1027/(1+exp((r-3.75)/0.52));
    }
    if(nucl == 3){
        dens = 0.1095/(1+exp((r-4.781)/0.52));
    }
    if(nucl == 4){
        dens = 0.1156/(1+exp((r-7.423)/0.52));
    }
    if(nucl == 5){
        dens = 0.1066/(1+exp((r-4.274)/0.52));
    }
    if(nucl == 6){
        dens = 0.1081/(1+exp((r-4.543)/0.52));
    }
    if(nucl == 7){
        dens = 0.1090/(1+exp((r-4.724)/0.52));
    }
    if(nucl == 8){
    
    } 
    if(nucl == 9){
        
    } 
    if(nucl == 10){
        
    }
    
    return dens;
}

// calculate transparency in z direction
double getTransparency(int nucl, double x, double y, double z, double rmax, double sigmaNN){
    // (x,y,z) = initial cell position
    double transp= 0.0;
    double dStep = 0.01;
    double rPrime;
    do { // integral rho sigma dz
        rPrime=std::sqrt(x*x+y*y+z*z);
        double dens = GenDensity(nucl, rPrime);
        transp = transp + (sigmaNN*dens*dStep);
        z = z + dStep;
    } while (rPrime<rmax);
    transp = exp(-transp);
    return transp;
}

double ClebshGordan(int l1, int l2){
    return 1.0;
}

// Associated Laguerre polynomial L_n^alpha(x)
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

// Normalization constant N_{n_r l}
double normalization(int n_r, int l, double b) {
    double num = 2.0 * std::tgamma(n_r + 1);
    double den = std::pow(b, 3) * std::tgamma(n_r + l + 1.5);
    return std::sqrt(num / den);
}

// Radial wavefunction R_{n_r l}(r)
double R_radial(int n_r, int l, double b, double r) {
    double rho = r / b;
    double rho2 = rho * rho;
    double alpha = l + 0.5;
    double L = assoc_laguerre(n_r, alpha, rho2);
    double N = normalization(n_r, l, b);
    return N * std::pow(rho, l) * L * std::exp(-rho2 / 2.0);
}

double PairProbability(){
    return pairProb;
}

// Cross section ratios ---------

double P = 0.05;
double sigmap = 2.4;
double sigman = 1;

double getReep(double np, double pp, double nn, double Tnp_p, double Tpp1, double Tpp2, double Tnp_n, double Tnn1, double Tnn2) {
    return (Tnp_p*np+(Tpp1+Tpp2)*pp)*sigmap + (Tnp_n*np+(Tnn1+Tnn2)*nn)*sigman*P;
}

int main(int argc, char* argv[]) {
    
    if(argc<4){
        std::cerr << "Usage: " << argv[0]
                  << " Nucleus type, 1- 12C, 2- 27Al, 3-Fe, 4-Pb, 5-40Ca, 6 - 48Ca\n" << "Rsrc pairProb";
        return 1;
    }
    pairProb = std::stod(argv[3]);
    
    // Initialize arrays
    for(int i=0;i<4;i++){
        for(int j=0;j<8;j++){
            for(int k=0;k<4;k++){
                for(int l=0;l<8;l++){
                    np[i][j][k][l]=0;
                    pp[i][j][k][l]=0;
                    nn[i][j][k][l]=0;
                }
            }
        }
    }
    
    const double LimitAngular = 100.;
    int nuclType = std::stoi(argv[1]);
    
    // Parameters
    const double Rsrc = std::stod(argv[2]);
    const double Rsrc_squared = Rsrc * Rsrc; // Precompute for faster comparison
    double r_max = 10; // fm
    const double spacing = 0.3; // fm
    const double sigmaNN = 3.5; //fm^2
    double A = 208.0;  // mass number
    
    // Nuclear configuration setup (same as original)
    if(nuclType==0)
        A = 2.;
    if(nuclType==1)
        A = 12.;
    if(nuclType==2)
        A = 27.;
    if(nuclType==3)
        A = 56.;
    if(nuclType==4)
        A = 208.;
    if(nuclType==5)
        A = 40.;
    if(nuclType==6)
        A = 48.;
    if(nuclType==7)
        A = 54.;
    if(nuclType==8)
        A = 9.;
    if(nuclType==9)
        A = 63.;
    if(nuclType==10)
        A = 197.;

    // Oscillator length b (ħω ≈ 41 A^(-1/3))
    const double hbarc = 197.3269804; // MeV·fm
    const double m_n   = 939.0;       // MeV
    double hbar_omega  = 41.0 * std::pow(A, -1.0/3.0);
    double b = hbarc / std::sqrt(m_n * hbar_omega);

    double hbar_omega_p;
    double hbar_omega_n;
    double b_p;
    double b_n;

    std::vector<ShellState> protons;
    std::vector<ShellState> neutrons;
    
    // Complete nucleus configuration
    if(nuclType==0){
        protons = {{0,0}};
        neutrons ={{0,0}};
        
        hbar_omega_p = 41.0 * std::pow(1, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(1, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    
    if(nuclType==1){ // 12Ca
        protons = {{0,0}, {0,0}, {0,1}, {0,1},{0,1},{0,1}};
        neutrons = protons;
        
        hbar_omega_p = 41.0 * std::pow(6, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(6, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    
    if(nuclType==2){ // 27Al
        protons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2}};
        neutrons = {{0,0},{0,0},{0,1},{0,1},{0,1},{0,1},{0,1},{0,1},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2}};
        
        hbar_omega_p = 41.0 * std::pow(13, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(14, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    
    if(nuclType==3){ // 56Fe
        protons ={        {0,0},
         {0,0},        {0,1},        {0,1},        {0,1},
         {0,1},        {0,1},        {0,1},        {0,2},
         {0,2},        {0,2},        {0,2},        {0,2},
         {0,2},        {0,2},        {0,2},        {0,2},
         {0,2},        {1,0},        {1,0},        {0,3},
         {0,3},        {0,3},        {0,3},        {0,3},
         {0,3} };
         
         neutrons ={        {0,0},
         {0,0},        {0,1},        {0,1},        {0,1},
         {0,1},        {0,1},        {0,1},        {0,2},
         {0,2},        {0,2},        {0,2},        {0,2},
         {0,2},        {0,2},        {0,2},        {0,2},
         {0,2},        {1,0},        {1,0},        {0,3},
         {0,3},        {0,3},        {0,3},        {0,3},
         {0,3},        {0,3},        {0,3},        {1,1},
         {1,1}};
        
        hbar_omega_p = 41.0 * std::pow(26, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(30, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    
    if(nuclType==4){
        // Build proton and neutron vectors for Pb-208
         // Define shells for protons and neutrons in 208Pb
         protons ={{0,0},{0,0},   // 1s1/2
         {0,1},{0,1},{0,1},{0,1},   // 1p3/2
         {0,1},{0,1},   // 1p1/2
         {0,2},{0,2},{0,2},{0,2},{0,2},{0,2},   // 1d5/2
         {1,0},{1,0},   // 2s1/2
         {0,2},{0,2},{0,2},{0,2},   // 1d3/2
         {0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},   // 1f7/2
         {1,1},{1,1},{1,1},{1,1},   // 2p3/2
         {0,3},{0,3},{0,3},{0,3},{0,3},{0,3},   // 1f5/2
         {1,1},{1,1},   // 2p1/2
         {0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},  // 1g9/2
         {1,2},{1,2},{1,2},{1,2},{1,2},{1,2},   // 2d5/2
         {0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},   // 1g7/2
         {2,0},{2,0},   // 3s1/2
         {1,2},{1,2},{1,2},{1,2},   // 2d3/2
         {0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5}};   // 1h11/2
         
         neutrons ={{0,0},{0,0},   // 1s1/2
         {0,1},{0,1},{0,1},{0,1},   // 1p3/2
         {0,1},{0,1},   // 1p1/2
         {0,2},{0,2},{0,2},{0,2},{0,2},{0,2},   // 1d5/2
         {1,0},{1,0},   // 2s1/2
         {0,2},{0,2},{0,2},{0,2},   // 1d3/2
         {0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},   // 1f7/2
         {1,1},{1,1},{1,1},{1,1},   // 2p3/2
         {0,3},{0,3},{0,3},{0,3},{0,3},{0,3},   // 1f5/2
         {1,1},{1,1},   // 2p1/2
         {0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},  // 1g9/2
         {1,2},{1,2},{1,2},{1,2},{1,2},{1,2},   // 2d5/2
         {0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},{0,4},   // 1g7/2
         {2,0},{2,0},   // 3s1/2
         {1,2},{1,2},{1,2},{1,2},   // 2d3/2
         {0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},{0,5},
         {1,3},{1,3},{1,3},{1,3},{1,3},{1,3},{1,3},{1,3},  // 2f7/2
         {2,1},{2,1},{2,1},{2,1},  // 3p3/2
         {1,3},{1,3},{1,3},{1,3},{1,3},{1,3},  // 2f5/2
         {2,1},{2,1},  // 3p1/2
         {0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6},{0,6}, // 1i13/2
         {2,4},{2,4},{2,4},{2,4},{2,4},{2,4},{2,4},{2,4},{2,4},{2,4}  // 2g9/2
         };
        
        hbar_omega_p = 41.0 * std::pow(82, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(126, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    if(nuclType==5){
        // Build proton and neutron vectors for Ca-40
        protons ={
            {0,0},{0,0},  // 1s1/2
            {0,1},{0,1},{0,1},{0,1},{0,1},{0,1},  // 1p3/2 + 1p1/2
            {0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},  // 1d5/2 + 1d3/2
            {1,0},{1,0}
        };
        
        neutrons ={
            {0,0},{0,0},  // 1s1/2
            {0,1},{0,1},{0,1},{0,1},{0,1},{0,1},  // 1p3/2 + 1p1/2
            {0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},  // 1d5/2 + 1d3/2
            {1,0},{1,0}
        };
        
        hbar_omega_p = 41.0 * std::pow(20, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(20, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    
    if(nuclType==6){
        // Build proton and neutron vectors for Ca-48
        protons ={
            {0,0},{0,0},  // 1s1/2
            {0,1},{0,1},{0,1},{0,1},{0,1},{0,1},  // 1p3/2 + 1p1/2
            {0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},  // 1d5/2 + 1d3/2
            {1,0},{1,0}
        };
        
        neutrons ={
            {0,0},{0,0},  // 1s1/2
            {0,1},{0,1},{0,1},{0,1},{0,1},{0,1},  // 1p3/2 + 1p1/2
            {0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},  // 1d5/2 + 1d3/2
            {1,0},{1,0},
            {0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3}
        };
        
        hbar_omega_p = 41.0 * std::pow(20, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(28, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    
    if(nuclType==7){
        // Build proton and neutron vectors for Fe-54
        protons ={
            {0,0},{0,0},  // 1s1/2
            {0,1},{0,1},{0,1},{0,1},{0,1},{0,1},  // 1p3/2 + 1p1/2
            {0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},  // 1d5/2 + 1d3/2
            {1,0},{1,0},  // 2s1/2
            {0,3},{0,3},{0,3},{0,3},{0,3},{0,3}  // 1f7/2 (partially filled)
        };
        
        neutrons ={
            {0,0},{0,0},
            {0,1},{0,1},{0,1},{0,1},{0,1},{0,1},
            {0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},{0,2},
            {1,0},{1,0},  // 2s1/2
            {0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3},{0,3}
        };
        
        hbar_omega_p = 41.0 * std::pow(26, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(28, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    if(nuclType==8){ // 9Be (Z=4, N=5)
        // Protons: 1s1/2 (2) + 1p3/2 (2) = 4 protons
        protons = {{0,0}, {0,0}, {0,1}, {0,1}};
        // Neutrons: 1s1/2 (2) + 1p3/2 (3) = 5 neutrons
        neutrons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}};
        
        hbar_omega_p = 41.0 * std::pow(4, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(5, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }

    if(nuclType==9){ // 63Cu (Z=29, N=34)
        // Protons: 1s1/2(2) + 1p3/2(4) + 1p1/2(2) + 1d5/2(6) + 2s1/2(2) + 1d3/2(4) + 1f7/2(8) + 2p3/2(1)
        protons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1},
                   {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {1,0}, {1,0},
                   {0,2}, {0,2}, {0,2}, {0,2}, {0,3}, {0,3}, {0,3}, {0,3},
                   {0,3}, {0,3}, {0,3}, {0,3}, {1,1}};
        
        // Neutrons: 1s1/2(2) + 1p3/2(4) + 1p1/2(2) + 1d5/2(6) + 2s1/2(2) + 1d3/2(4) + 1f7/2(8) + 2p3/2(4) + 1f5/2(2)
        neutrons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1},
                    {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {1,0}, {1,0},
                    {0,2}, {0,2}, {0,2}, {0,2}, {0,3}, {0,3}, {0,3}, {0,3},
                    {0,3}, {0,3}, {0,3}, {0,3}, {1,1}, {1,1}, {1,1}, {1,1},
                    {0,3}, {0,3}};
        
        hbar_omega_p = 41.0 * std::pow(29, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(34, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }
    
    if(nuclType==10){ // 197Au (Z=79, N=118)
        // Protons: Complete shells up to 2d3/2(1) for total of 79
        protons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1},
                   {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {1,0}, {1,0},
                   {0,2}, {0,2}, {0,2}, {0,2}, {0,3}, {0,3}, {0,3}, {0,3},
                   {0,3}, {0,3}, {0,3}, {0,3}, {1,1}, {1,1}, {1,1}, {1,1},
                   {0,3}, {0,3}, {1,1}, {1,1}, {2,0}, {2,0}, {1,2}, {1,2},
                   {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2},
                   {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4},
                   {0,4}, {0,4}, {0,4}, {0,4}, {1,3}, {1,3}, {1,3}, {1,3},
                   {1,3}, {1,3}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1},
                   {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1},
                   {2,2}, {2,2}, {2,2}};
        
        // Neutrons: Complete shells up to 1i13/2(14) for total of 118
        neutrons = {{0,0}, {0,0}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1},
                    {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {0,2}, {1,0}, {1,0},
                    {0,2}, {0,2}, {0,2}, {0,2}, {0,3}, {0,3}, {0,3}, {0,3},
                    {0,3}, {0,3}, {0,3}, {0,3}, {1,1}, {1,1}, {1,1}, {1,1},
                    {0,3}, {0,3}, {1,1}, {1,1}, {2,0}, {2,0}, {1,2}, {1,2},
                    {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2}, {1,2},
                    {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4}, {0,4},
                    {0,4}, {0,4}, {0,4}, {0,4}, {1,3}, {1,3}, {1,3}, {1,3},
                    {1,3}, {1,3}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1},
                    {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1}, {2,1},
                    {2,2}, {2,2}, {2,2}, {2,2}, {2,2}, {2,2}, {1,4}, {1,4},
                    {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4},
                    {2,3}, {2,3}, {2,3}, {2,3}, {3,0}, {3,0}, {2,2}, {2,2},
                    {2,2}, {2,2}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4}, {1,4},
                    {1,4}, {1,4}, {1,4}, {1,4}, {2,3}, {2,3}, {2,3}, {2,3},
                    {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5},
                    {1,5}, {1,5}, {1,5}, {1,5}, {1,5}, {1,5}};
        
        hbar_omega_p = 41.0 * std::pow(79, -1.0/3.0);
        hbar_omega_n = 41.0 * std::pow(118, -1.0/3.0);
        b_p =  hbarc / std::sqrt(m_n * hbar_omega_p);
        b_n =  hbarc / std::sqrt(m_n * hbar_omega_n);
    }

    // Pre-calculate centers
    std::vector<std::array<double,3>> centers;
    const int n_steps_1d = static_cast<int>(std::floor(2 * r_max / spacing)) + 1;
    centers.reserve(n_steps_1d * n_steps_1d * n_steps_1d);
    
    for (int i = 0; i < n_steps_1d; ++i) {
        double x = -r_max + i * spacing;
        for (int j = 0; j < n_steps_1d; ++j) {
            double y = -r_max + j * spacing;
            for (int k = 0; k < n_steps_1d; ++k) {
                double z = -r_max + k * spacing;
                centers.push_back({x, y, z});
            }
        }
    }

    // After generating the `centers` vector, precompute useful per-center quantities
    double cellVol = std::pow(spacing,3);

    // --- Optimization: precompute radius and transparency for every center ---
    std::vector<double> center_r(centers.size());
    std::vector<double> center_transparency(centers.size());

    #pragma omp parallel for schedule(dynamic, 1000)
    for (size_t idx = 0; idx < centers.size(); ++idx) {
        const auto &c = centers[idx];
        double r = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
        center_r[idx]           = r;
        center_transparency[idx] = getTransparency(nuclType, c[0], c[1], c[2], r_max, sigmaNN);
    }
    // ------------------------------------------------------------------------

    std::vector<Nucleon> nucleons;
    nucleons.reserve(protons.size() * centers.size() + neutrons.size() * centers.size());
    
    std::cout << "Creating nucleon probability distribution..." << std::endl;
    
    // Probability threshold
    const double accuracy = 1e-3;
    double total_volume = std::pow(2.0 * r_max, 3);
    double N_nucleons = A * total_volume/std::pow(spacing, 3);
    double N_pairs = N_nucleons*(N_nucleons-1) / 2;
    const double prob_threshold = std::sqrt( accuracy / N_pairs );
    cout << "Probability threshold: " << prob_threshold << endl;
    
    // Parallelize nucleon creation
    std::vector<std::vector<Nucleon>> thread_nucleons(omp_get_max_threads());
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        
        // Add protons at all possible centers
        #pragma omp for schedule(dynamic, 1)
        for (size_t p_idx = 0; p_idx < protons.size(); ++p_idx) {
            auto& p = protons[p_idx];
            for (size_t c_idx = 0; c_idx < centers.size(); ++c_idx) {
                auto& c = centers[c_idx];
                double r = center_r[c_idx];
                double Rp = R_radial(p.n_r, p.l, b_p, r);
                double prob = (Rp*Rp/(4*M_PI)) * cellVol;
                
                if (prob > prob_threshold) {
                    double transp = center_transparency[c_idx];
                    thread_nucleons[thread_id].push_back({0, (int)p_idx, c_idx, prob, transp});
                }
            }
        }
        
        // Add neutrons at all possible centers
        #pragma omp for schedule(dynamic, 1)
        for (size_t n_idx = 0; n_idx < neutrons.size(); ++n_idx) {
            auto& n = neutrons[n_idx];
            for (size_t c_idx = 0; c_idx < centers.size(); ++c_idx) {
                auto& c = centers[c_idx];
                double r = center_r[c_idx];
                double Rn = R_radial(n.n_r, n.l, b_n, r);
                double prob = (Rn*Rn/(4*M_PI)) * cellVol;
                
                if (prob > prob_threshold) {
                    double transp = center_transparency[c_idx];
                    thread_nucleons[thread_id].push_back({1, (int)n_idx, c_idx, prob, transp});
                }
            }
        }
    }
    
    // Combine results from all threads
    for (auto& tn : thread_nucleons) {
        nucleons.insert(nucleons.end(), tn.begin(), tn.end());
    }
    
    std::cout << "Created " << nucleons.size() << " nucleon-position combinations" << std::endl;

    // Thread-local storage for results
    std::vector<double> thread_P_np(omp_get_max_threads(), 0.0);
    std::vector<double> thread_P_pp(omp_get_max_threads(), 0.0);
    std::vector<double> thread_P_nn(omp_get_max_threads(), 0.0);
    std::vector<double> thread_Reep(omp_get_max_threads(), 0.0);
    
    // Thread-local arrays for nuclear pair accumulation (flattened for performance)
    const int n_r_max = 4;
    const int l_max = 8;
    const size_t shell_state_size = (size_t)n_r_max * l_max * n_r_max * l_max;

    std::vector<double> thread_np(omp_get_max_threads() * shell_state_size, 0.0);
    std::vector<double> thread_pp(omp_get_max_threads() * shell_state_size, 0.0);
    std::vector<double> thread_nn(omp_get_max_threads() * shell_state_size, 0.0);
    
    std::cout << "Computing pair probabilities..." << std::endl;
    
    // --- Spatial Partitioning Optimization ---
    // Group nucleons by their center index using a flattened structure for better data locality
    std::vector<std::vector<size_t>> temp_nucleons_in_cell(centers.size());
    for (size_t i = 0; i < nucleons.size(); ++i) {
        temp_nucleons_in_cell[nucleons[i].center_index].push_back(i);
    }

    std::vector<size_t> cell_offsets(centers.size() + 1, 0);
    size_t total_nucleons_in_cells = 0;
    for(size_t i = 0; i < centers.size(); ++i) {
        total_nucleons_in_cells += temp_nucleons_in_cell[i].size();
        cell_offsets[i+1] = total_nucleons_in_cells;
    }

    std::vector<size_t> flat_nucleons_in_cell;
    flat_nucleons_in_cell.reserve(total_nucleons_in_cells);
    for(size_t i = 0; i < centers.size(); ++i) {
        flat_nucleons_in_cell.insert(flat_nucleons_in_cell.end(), temp_nucleons_in_cell[i].begin(), temp_nucleons_in_cell[i].end());
    }
    temp_nucleons_in_cell.clear(); // Free memory

    const int search_radius_cells = static_cast<int>(ceil(Rsrc / spacing));
    
    // Parallelize the main computation loop
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        
        auto process_pair = [&](const Nucleon& n1, const Nucleon& n2) {
            // Initialize transparencies
            double Tnp_p = 0.0, Tnp_n = 0.0;
            double Tpp1 = 0.0, Tpp2 = 0.0;
            double Tnn1 = 0.0, Tnn2 = 0.0;
            double np_prob = 0.0, pp_prob = 0.0, nn_prob = 0.0;
            
            // Check if it's the same nucleon
            if (n1.type == n2.type && n1.shell_index == n2.shell_index) return;
            
            // Calculate pair probability
            double pair_prob = n1.probability * n2.probability;
            
            // n-p pairs
            if ((n1.type == 0 && n2.type == 1) || (n1.type == 1 && n2.type == 0)) {
                int p_idx = (n1.type == 0) ? n1.shell_index : n2.shell_index;
                int n_idx = (n1.type == 1) ? n1.shell_index : n2.shell_index;
                auto& p = protons[p_idx];
                auto& n = neutrons[n_idx];
                
                if (n1.type == 0) {
                    Tnp_p = n1.transparency;
                    Tnp_n = n2.transparency;
                } else {
                    Tnp_p = n2.transparency;
                    Tnp_n = n1.transparency;
                }
                
                if (std::abs(p.l - n.l) <= LimitAngular) {
                    np_prob = pair_prob;
                    thread_P_np[thread_id] += np_prob;
                    size_t flat_idx = thread_id * shell_state_size + (size_t)n.n_r * (l_max * n_r_max * l_max) + (size_t)n.l * (n_r_max * l_max) + (size_t)p.n_r * l_max + (size_t)p.l;
                    thread_np[flat_idx] += np_prob;
                }
            }
            // p-p pairs
            else if (n1.type == 0 && n2.type == 0) {
                auto& p1 = protons[n1.shell_index];
                auto& p2 = protons[n2.shell_index];
                
                Tpp1 = n1.transparency;
                Tpp2 = n2.transparency;
                
                if (std::abs(p1.l - p2.l) <= LimitAngular) {
                    pp_prob = pair_prob * pairProb;
                    thread_P_pp[thread_id] += pp_prob;
                    size_t flat_idx = thread_id * shell_state_size + (size_t)p1.n_r * (l_max * n_r_max * l_max) + (size_t)p1.l * (n_r_max * l_max) + (size_t)p2.n_r * l_max + (size_t)p2.l;
                    thread_pp[flat_idx] += pp_prob;
                }
            }
            // n-n pairs
            else if (n1.type == 1 && n2.type == 1) {
                auto& n1_state = neutrons[n1.shell_index];
                auto& n2_state = neutrons[n2.shell_index];
                
                Tnn1 = n1.transparency;
                Tnn2 = n2.transparency;
                
                if (std::abs(n1_state.l - n2_state.l) <= LimitAngular) {
                    nn_prob = pair_prob * pairProb;
                    thread_P_nn[thread_id] += nn_prob;
                    size_t flat_idx = thread_id * shell_state_size + (size_t)n1_state.n_r * (l_max * n_r_max * l_max) + (size_t)n1_state.l * (n_r_max * l_max) + (size_t)n2_state.n_r * l_max + (size_t)n2_state.l;
                    thread_nn[flat_idx] += nn_prob;
                }
            }
            
            thread_Reep[thread_id] += getReep(np_prob, pp_prob, nn_prob, Tnp_p, Tpp1, Tpp2, Tnp_n, Tnn1, Tnn2);
        };

        #pragma omp for schedule(dynamic, 1)
        for (int ix1 = 0; ix1 < n_steps_1d; ++ix1) {
            for (int iy1 = 0; iy1 < n_steps_1d; ++iy1) {
                for (int iz1 = 0; iz1 < n_steps_1d; ++iz1) {
                    size_t c1_idx = (size_t)ix1 * n_steps_1d * n_steps_1d + (size_t)iy1 * n_steps_1d + iz1;

                    size_t start1 = cell_offsets[c1_idx];
                    size_t end1 = cell_offsets[c1_idx + 1];
                    if (start1 == end1) continue;

                    auto& c1 = centers[c1_idx];

                    // Pairs within the same cell
                    for (size_t i = start1; i < end1; ++i) {
                        for (size_t j = i + 1; j < end1; ++j) {
                             process_pair(nucleons[flat_nucleons_in_cell[i]], nucleons[flat_nucleons_in_cell[j]]);
                        }
                    }

                    // Pairs with neighboring cells
                    for (int ix2 = ix1 - search_radius_cells; ix2 <= ix1 + search_radius_cells; ++ix2) {
                        for (int iy2 = iy1 - search_radius_cells; iy2 <= iy1 + search_radius_cells; ++iy2) {
                            for (int iz2 = iz1 - search_radius_cells; iz2 <= iz1 + search_radius_cells; ++iz2) {
                                if (ix2 < 0 || ix2 >= n_steps_1d || iy2 < 0 || iy2 >= n_steps_1d || iz2 < 0 || iz2 >= n_steps_1d) continue;

                                size_t c2_idx = (size_t)ix2 * n_steps_1d * n_steps_1d + (size_t)iy2 * n_steps_1d + iz2;
                                if (c2_idx <= c1_idx) continue;

                                size_t start2 = cell_offsets[c2_idx];
                                size_t end2 = cell_offsets[c2_idx + 1];
                                if (start2 == end2) continue;

                                auto& c2 = centers[c2_idx];
                                double dx = c1[0] - c2[0];
                                double dy = c1[1] - c2[1];
                                double dz = c1[2] - c2[2];
                                double dist_squared = dx * dx + dy * dy + dz * dz;

                                if (dist_squared >= Rsrc_squared) continue;
                                
                                for (size_t n1_flat_idx = start1; n1_flat_idx < end1; ++n1_flat_idx) {
                                    for (size_t n2_flat_idx = start2; n2_flat_idx < end2; ++n2_flat_idx) {
                                        process_pair(nucleons[flat_nucleons_in_cell[n1_flat_idx]], nucleons[flat_nucleons_in_cell[n2_flat_idx]]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Combine results from all threads
    double P_np = 0.0, P_pp = 0.0, P_nn = 0.0, Reep = 0.0;
    for (int t = 0; t < omp_get_max_threads(); ++t) {
        P_np += thread_P_np[t];
        P_pp += thread_P_pp[t];
        P_nn += thread_P_nn[t];
        Reep += thread_Reep[t];
        
        // Combine arrays
        size_t thread_offset = t * shell_state_size;
        for (int i = 0; i < n_r_max; ++i) {
            for (int j = 0; j < l_max; ++j) {
                for (int k = 0; k < n_r_max; ++k) {
                    for (int l = 0; l < l_max; ++l) {
                        size_t flat_idx = (size_t)i * (l_max * n_r_max * l_max) + (size_t)j * (n_r_max * l_max) + (size_t)k * l_max + l;
                        np[i][j][k][l] += thread_np[thread_offset + flat_idx];
                        pp[i][j][k][l] += thread_pp[thread_offset + flat_idx];
                        nn[i][j][k][l] += thread_nn[thread_offset + flat_idx];
                    }
                }
            }
        }
    }
    
    // Calculate total probabilities
    double Pp = 0.0, Pn = 0.0;
    for (auto& n : nucleons) {
        if (n.type == 0) Pp += n.probability;
        else Pn += n.probability;
    }

    // Output results
    std::cout << "\nTotal close-proximity probabilities (small spheres):\n";
    std::cout << "  #np: " << P_np << "\n";
    std::cout << "  #pp: " << P_pp << "\n";
    std::cout << "  #nn: " << P_nn << "\n";
    std::cout << "  R(e,e'p): " << Reep << "\n";
    
    cout << "rmax: " << r_max << "\n";
    cout << " SRC fraction = " << (P_np+P_pp+P_nn)/A << endl;
    cout << " pp/np = " << P_pp/P_np << "\t nn/np = " << P_nn/P_np << endl;
    cout << " Pp = " << Pp << " Pn " << Pn << endl;
    cout << "spacing: " << spacing << "\n";
    cout << " Rsrc = " << Rsrc << endl;
    cout << " pairProb = " << pairProb << endl;
    
    return 0;
}
