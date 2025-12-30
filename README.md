# Stochastic Modeling of Stock Prices based on Geometric Brownian Motion (GBM)

## Project Overview
This repository contains the MATLAB source code for the "Applied Stochastic Processes" course project. The project simulates stock price trajectories using Geometric Brownian Motion (GBM) and validates the properties of Itô calculus.

**Course:** Applied Stochastic Processes  
**Language:** MATLAB  
**Year:** 2025

## Code Structure (`GBM_Simulation.m`)
The script includes the following modules:
1.  **Parameter Settings:** Initialization of drift ($\mu$), volatility ($\sigma$), and time steps.
2.  **Path Simulation:** Monte Carlo simulation of stock price paths using the analytical solution of GBM.
3.  **Sensitivity Analysis:** Visualizing the impact of volatility parameters on asset prices.
4.  **Quadratic Variation:** Numerical verification of the convergence to $\sigma^2 t$.
5.  **Empirical Analysis:** Parameter estimation (MLE) and simulation based on real market data logic.

## How to Run
1. Ensure MATLAB is installed.
2. Download `main.m`.
3. Run the script directly. It will generate 4 figures automatically:
    - Figure 1: Simulated Paths
    - Figure 2: Sensitivity Analysis ($\sigma=0.1$ vs $\sigma=0.6$)
    - Figure 3: Convergence of Quadratic Variation
    - Figure 4: Empirical Analysis vs GBM Forecast

## Author
[李尚嵋]
[125034910008]
