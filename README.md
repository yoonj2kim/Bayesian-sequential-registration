# Bayesian Sequential Registration for Functional Data
MATLAB code for implementation of [Bayesian sequential registration for functional data](https://arxiv.org/abs/2203.12005) by Yoonji Kim, Oksana Chkrebtii, and Sebastian Kurtek. To publish work based on this code, please cite: Yoonji Kim, Oksana Chkrebtii, and Sebastian Kurtek, “Sequential Bayesian Registration for Functional Data,” arXiv:2203.12005 [stat.ME], Mar. 2022.

## Set up
Before running the code, please download [the fdasrvf package](https://github.com/jdtuck/fdasrvf_MATLAB) and [the constrained piecewise linear least square fit package](https://www.mathworks.com/matlabcentral/fileexchange/40913-piecewise-linear-least-square-fit).

## Code description:

- ['simulate_example.m'](simulate_example.m): This code generates data and results as demonstrated in Section 5.2.
- ['mcmc_regist.m'](mcmc_regist.m): This function initializes posterior samples for function registration using a MCMC sampler.
- ['smc_regist.m'](smc_regist.m): Given existing posterior particles, this function allows us to update the posterior particles when a new function arrives.

## Datasets:

- 'simulated_example_1_data.mat': Simulated data and estimation results presented in Section 5.2.
- 'simulated_example_2_data.mat': Simulated data and estimation results presented in Section 5.3.
- 'real_data_analysis_1_data.mat': Real data analysis on drought intensity of Kaweah river in California (Section 5.4).
- 'real_data_analysis_2_elnino_data.mat', 'real_data_analysis_2_neutral_data.mat', 'real_data_analysis_2_lanina_data.mat': Real data analysis on sea surface salinity in Pacific Ocean and El Niño-Southern Oscillation (Section 5.5).
- 'real_data_analysis_3_data.mat': Real data analysis on sequence of repeated patterns in electrocardiogram signals (Section 5.6).

# Reference

Guido Albertin (2022). Piecewise linear least square fit (https://www.mathworks.com/matlabcentral/fileexchange/40913-piecewise-linear-least-square-fit), MATLAB Central File Exchange. Retrieved October 18, 2022.

Yoonji Kim, Oksana Chkrebtii, and Sebastian Kurtek, “Sequential Bayesian Registration for Functional Data,” arXiv:2203.12005 [stat.ME], Mar. 2022.

tetonedge (2022). fdasrvf (https://github.com/jdtuck/fdasrvf_MATLAB), GitHub. Retrieved October 18, 2022.
