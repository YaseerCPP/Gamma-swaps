#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>

// Function to simulate asset price paths using the Black-Scholes model
std::vector<std::vector<double>> simulatePricePaths(double S0, double mu, double sigma, double T, int numSteps, int numPaths) {
    std::vector<std::vector<double>> pricePaths(numPaths, std::vector<double>(numSteps + 1, S0));
    double dt = T / numSteps;
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<> distribution(0.0, 1.0);

    for (int i = 0; i < numPaths; ++i) {
        for (int j = 1; j <= numSteps; ++j) {
            double dW = distribution(generator) * sqrt(dt);
            pricePaths[i][j] = pricePaths[i][j-1] * exp((mu - 0.5 * sigma * sigma) * dt + sigma * dW);
        }
    }
    return pricePaths;
}

// Function to calculate the gamma of a European call option
double calculateGamma(double S, double K, double T, double r, double sigma) {
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double gamma = exp(-0.5 * d1 * d1) / (S * sigma * sqrt(2 * M_PI * T));
    return gamma;
}

// Function to calculate the realized gamma over multiple price paths
double calculateRealizedGamma(const std::vector<std::vector<double>>& pricePaths, double K, double T, double r, double sigma) {
    double totalGamma = 0.0;
    int count = 0;
    int numSteps = pricePaths[0].size() - 1;
    double dt = T / numSteps;

    for (const auto& path : pricePaths) {
        for (int j = 1; j <= numSteps; ++j) {
            double gamma = calculateGamma(path[j], K, T - j * dt, r, sigma);
            totalGamma += gamma;
            ++count;
        }
    }
    return totalGamma / count;
}

// Main function to price the gamma swap using Monte Carlo simulation
int main() {
    // Parameters
    double S0 = 100.0;    // Initial asset price
    double K = 100.0;     // Strike price
    double T = 1.0;       // Time to maturity in years
    double r = 0.05;      // Risk-free interest rate
    double sigma = 0.2;   // Volatility
    double mu = 0.1;      // Drift
    int numSteps = 252;   // Number of time steps (e.g., daily steps)
    int numPaths = 10000; // Number of Monte Carlo paths

    // Simulate price paths
    auto pricePaths = simulatePricePaths(S0, mu, sigma, T, numSteps, numPaths);

    // Calculate the realized gamma
    double realizedGamma = calculateRealizedGamma(pricePaths, K, T, r, sigma);

    // Example strike gamma (expected gamma)
    double strikeGamma = 0.01; // Example value

    // Calculate the payoff
    double notional = 1000000; // Example notional amount
    double payoff = (realizedGamma - strikeGamma) * notional;

    // Output the result
    std::cout << "Realized Gamma: " << realizedGamma << std::endl;
    std::cout << "Payoff: " << payoff << std::endl;

    return 0;
}
