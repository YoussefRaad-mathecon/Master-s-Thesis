# Master's Thesis: The Black-Scholes Option Pricing Model - A Markov-Switching Extension

## Overview
This Master's thesis explores extensions to the Black-Scholes option pricing model as a frequent critique of said model is it lack of adaptivity to economical turbulent environments. The extensions are via continuous and discrete state-space models often called Markov regime-switching models. The analysis is done using S\&P 500 price data and total return data.

**Practical Information**
- **Author**: Youssef Raad
- **Supervisor**: Rolf Poulsen
- **University**: University of Copenhagen
- **Department**: Mathmematical Science
- **Date**: 22<sup>nd</sup> December 2025

## Abstract
The Black-Scholes model (BSM) has long been a cornerstone of financial theory; however, its assumption of constant drift and volatility fails to capture the time-varying nature of asset returns, such as volatility clustering and abrupt regime shifts. This thesis investigates whether extending the BSM to include dynamic parameter evolution improves performance. Two extensions are proposed and implemented: a Black-Scholes Hidden Markov Model (BS-HMM) and a Black-Scholes Continuous State-Space Model (BS-SSM).

Using daily S\&P 500 data from 1927 to 2025, the models are calibrated via Maximum Likelihood Estimation. In-sample analysis reveals that the extended models, particularly a 4-state BS-HMM and a factor-loaded continuous state-space model ($\text{BS-SSM}_{\beta}$), provide a superior fit to historical data compared to the static BSM. These models successfully identify distinct market phases, distinguishing between tranquil bull markets and high-volatility crisis regimes, such as the 1929 Crash and the 2008 Financial Crisis.

Despite the richer descriptive power and improved in-sample fit, out-of-sample evaluation on a hold-out period (2020--2025) indicates that the regime-switching extensions offer negligible gains in one-step-ahead point forecasting accuracy (MSE and RMSE) relative to the constant-parameter BSM. While the extended models produce more realistic, horizon-dependent forecast densities, the findings suggest that the added complexity of latent state inference does not translate into superior short-term predictive power for point forecasts.

## Repository Structure
This section illustrates the structure of code (in `Python` and `R`) and dataframes (in `Excel` and `R`). The total amount of files is 53.
```
├── README.md                          
├── CodePython/                              
│   ├── DividendsPython.py
│   └── SP500DataGeneratorPython.py
├── CodeR/                              
│   ├── AICBIC.R
│   ├── AICBICBSHMM.R
│   ├── AICBICBSSSM.R
│   ├── BSHMMPredictionInterval.R
│   ├── BSMPredictionInterval.R
│   ├── BSSSMBetaPredictionInterval.R
│   ├── DataR.R
│   ├── DividendsR.R
│   ├── DividendStatewise.R
│   ├── FittingBSHMM.R
│   ├── FittingBSM.R
│   ├── FittingBSMResidualsStandardErrors.R
│   ├── FittingBSSSM.R
│   ├── FittingBSSSMBeta.R
│   ├── FittingBSSSMRobustness.R
│   ├── FittingBSSSMRobustnessBeta.R
│   ├── NumberOfParameters.R
│   ├── PlottingBSM.R
│   ├── PlottingStateDependentDistributionsBSHMM.R
│   ├── PseudoResidualsBSHMM.R
│   ├── PseudoResidualsBSSSM.R
│   ├── PseudoResidualsBSSSMBeta.R
│   ├── SimulationStudyBS.R
│   ├── SimulationStudyBSHMM.R
│   ├── SimulationStudyBSSSM.R
│   ├── SimulationStudyBSSSMBeta.R
│   ├── StandardErrorsBSHMM.R
│   ├── StandardErrorsBSSSM.R
│   ├── StandardErrorsBSSSMBeta.R
│   ├── ViterbiBSSSMBeta.R
│   └── ViterbiHMM.R
├── DataFramesExcel/                                
│   ├── sp500_clean_with_divs.csv
│   ├── sp500_clean_with_divs_1927_2019.csv
│   ├── sp500_clean_with_divs_FULL.csv
│   ├── sp500_dividend_daily_MATCHED_to_gspc.csv
│   ├── sp500_dividend_daily_backfilled.csv
│   ├── sp500_yahoo_daily_full.csv
│   ├── statewise_dividend_yields_4and5state_with_SEs.csv          
│   └── statewise_dividend_yields_with_SEs.csv
├── DataFramesR/                         
│   ├── BSSSM_loglik_comparison.RData        
│   ├── fitted_params_BSHMM_2019.RData 
│   ├── fitted_params_BSSSM_2019.RData          
│   ├── fitted_params_BSSSM_2019_beta.RData      
│   ├── fitted_params_BSSSM_2019_test.RData
│   ├── fitted_params_BSSSM_grid_2019.RData
│   ├── fitted_params_BSSSM_grid_2019_beta.RData     
│   ├── fitted_params_BS_2019.RData
│   ├── sp500_clean_with_divs.RData
│   ├── sp500_clean_with_divs_1927_2019.RData
│   └── sp500_clean_with_divs_FULL.RData
```
