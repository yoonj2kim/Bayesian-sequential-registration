
% setup paths
addpath(genpath('DP'))
addpath(genpath('bspline_tools'))
addpath(genpath('basis'))

%% El Nino years

%% Set up
load('data/real_data_analysis_2_elnino_data.mat','f_mat','t','tG','M','N','M_gamma','nbasis','q_basis','prior_params')

q_mat = zeros(M,N);
for i = 1:N
    q_mat(:,i) = f_to_srvf(f_mat(:,i),t);
end

%% Initialization
% Derive posterior samples given functions 1 through
% N_func using MCMC.

N_func = 3;

% initial samples
c_ini = zeros(1,nbasis);
ds_ini = repmat(tG',[N_func 1]);
sigma_squared_ini = 0.1;

% set hyperparameters and tuning parameters
kappa_prop = 50*ones(N_func,1); % proposal concentration parameter for phase components
prior_params.kappa = 5; % prior concentration parameter for phase components
prior_params.c_mean = zeros(nbasis,1);
prior_params.c_var = 20*eye(nbasis); % prior mean and variance for basis coefficients of template function
prior_params.alpha = 4; prior_params.beta = 0.01; % hyperparameters for \sigma^2
stepvar_chain_ini = 0.05*cov(((q_basis'*q_basis)\q_basis'*q_mat)'); % initial proposal variance of 'c'
checknum_stop_fix = 1000; % the number of proposal variance adjustment
div_iter_by = 1000; % proposal variance adjustment every 'div_iter_by' number of MCMC iterations
target_rate = [0.2 0.55]; % target MH acceptance rate

qs = q_mat(:,1:N_func); % given SRVFs
total_mcmc = 10000; % total number of MCMC iterations

% run MCMC
[c_samples, d_samples, sigma_samples, log_post, kappa_prop, ...
    MCMC_time] = mcmc_regist(total_mcmc, c_ini, ds_ini, ...
    sigma_squared_ini, t, qs, q_basis, prior_params, kappa_prop, ...
    stepvar_chain_ini, target_rate, checknum_stop_fix,div_iter_by);

%% Sequential Bayesian Updating
% Assign equal weights to posterior samples obtained using MCMC
% Update weighted posterior samples through SMC when a new function arrives

J_particles = 10000; % number of posterior particles

% burn-in MCMC posterior samples and assign equal weights
c = c_samples((end-J_particles+1):end,:);
ds = d_samples((end-J_particles+1):end,:,:);
sigma_squared = sigma_samples((end-J_particles+1):end);
weights = 1/J_particles*ones(1,J_particles);

% proposal concentration parameter for phase components
kappa_prop = [kappa_prop(1:N_func);mean(kappa_prop(1:N_func))*ones(N-N_func,1)];

c_smc = cell(N,1);
d_smc = cell(N,1);
sigma_squared_smc = zeros(J_particles,N);
w_smc = zeros(J_particles,N);
ESS_smc = zeros(2,N);
logpost_smc = zeros(J_particles,N);
SMC_timerecord = zeros(N,1);

K_mh_c = 30;
K_mh_d = 30;

for N_func = (N_func+1):N

    qs = q_mat(:,1:N_func); % given SRVFs

    % run SMC algorithm
    [c_update, d_update, sigma_squared_update, w_update,...
        ESS, logpost_update, SMC_time] = smc_regist(c, ds, ...
        sigma_squared, weights, t, qs, q_basis, ...
        prior_params, kappa_prop, K_mh_c, K_mh_d);

    c_smc{N_func} = c_update;
    d_smc{N_func} = d_update;
    sigma_squared_smc(:,N_func) = sigma_squared_update;
    w_smc(:,N_func) = w_update;
    ESS_smc(:,N_func) = ESS;
    logpost_smc(:,N_func) = logpost_update;
    SMC_timerecord(N_func) = SMC_time;
    
    c = c_update;
    ds = d_update;
    sigma_squared = sigma_squared_update;
    weights = w_update;
end

%% Neutral years

%% Set up
load('real_data_analysis_2_neutral_data.mat','f_mat','t','tG','M','N','M_gamma','nbasis','q_basis','prior_params')

q_mat = zeros(M,N);
for i = 1:N
    q_mat(:,i) = f_to_srvf(f_mat(:,i),t);
end

%% Initialization
% Derive posterior samples given functions 1 through
% N_func using MCMC.

N_func = 3;

% initial samples
c_ini = zeros(1,nbasis);
ds_ini = repmat(tG',[N_func 1]);
sigma_squared_ini = 0.1;

% set hyperparameters and tuning parameters
kappa_prop = 50*ones(N_func,1); % proposal concentration parameter for phase components
prior_params.kappa = 5; % prior concentration parameter for phase components
prior_params.c_mean = zeros(nbasis,1);
prior_params.c_var = 20*eye(nbasis); % prior mean and variance for basis coefficients of template function
prior_params.alpha = 4; prior_params.beta = 0.01; % hyperparameters for \sigma^2
stepvar_chain_ini = 0.05*cov(((q_basis'*q_basis)\q_basis'*q_mat)'); % initial proposal variance of 'c'
checknum_stop_fix = 1000; % the number of proposal variance adjustment
div_iter_by = 1000; % proposal variance adjustment every 'div_iter_by' number of MCMC iterations
target_rate = [0.2 0.55]; % target MH acceptance rate

qs = q_mat(:,1:N_func); % given SRVFs
total_mcmc = 10000; % total number of MCMC iterations

% run MCMC
[c_samples, d_samples, sigma_samples, log_post, kappa_prop, ...
    MCMC_time] = mcmc_regist(total_mcmc, c_ini, ds_ini, ...
    sigma_squared_ini, t, qs, q_basis, prior_params, kappa_prop, ...
    stepvar_chain_ini, target_rate, checknum_stop_fix,div_iter_by);

%% Sequential Bayesian Updating
% Assign equal weights to posterior samples obtained using MCMC
% Update weighted posterior samples through SMC when a new function arrives

J_particles = 10000; % number of posterior particles

% burn-in MCMC posterior samples and assign equal weights
c = c_samples((end-J_particles+1):end,:);
ds = d_samples((end-J_particles+1):end,:,:);
sigma_squared = sigma_samples((end-J_particles+1):end);
weights = 1/J_particles*ones(1,J_particles);

% proposal concentration parameter for phase components
kappa_prop = [kappa_prop(1:N_func);mean(kappa_prop(1:N_func))*ones(N-N_func,1)];

c_smc = cell(N,1);
d_smc = cell(N,1);
sigma_squared_smc = zeros(J_particles,N);
w_smc = zeros(J_particles,N);
ESS_smc = zeros(2,N);
logpost_smc = zeros(J_particles,N);
SMC_timerecord = zeros(N,1);

K_mh_c = 30;
K_mh_d = 30;

for N_func = (N_func+1):N

    qs = q_mat(:,1:N_func); % given SRVFs

    % run SMC algorithm
    [c_update, d_update, sigma_squared_update, w_update,...
        ESS, logpost_update, SMC_time] = smc_regist(c, ds, ...
        sigma_squared, weights, t, qs, q_basis, ...
        prior_params, kappa_prop, K_mh_c, K_mh_d);

    c_smc{N_func} = c_update;
    d_smc{N_func} = d_update;
    sigma_squared_smc(:,N_func) = sigma_squared_update;
    w_smc(:,N_func) = w_update;
    ESS_smc(:,N_func) = ESS;
    logpost_smc(:,N_func) = logpost_update;
    SMC_timerecord(N_func) = SMC_time;
    
    c = c_update;
    ds = d_update;
    sigma_squared = sigma_squared_update;
    weights = w_update;
end


%% La Nina years

%% Set up
load('real_data_analysis_2_lanina_data.mat','f_mat','t','tG','M','N','M_gamma','nbasis','q_basis','prior_params')

q_mat = zeros(M,N);
for i = 1:N
    q_mat(:,i) = f_to_srvf(f_mat(:,i),t);
end

%% Initialization
% Derive posterior samples given functions 1 through
% N_func using MCMC.

N_func = 3;

% initial samples
c_ini = zeros(1,nbasis);
ds_ini = repmat(tG',[N_func 1]);
sigma_squared_ini = 0.1;

% set hyperparameters and tuning parameters
kappa_prop = 50*ones(N_func,1); % proposal concentration parameter for phase components
prior_params.kappa = 5; % prior concentration parameter for phase components
prior_params.c_mean = zeros(nbasis,1);
prior_params.c_var = 20*eye(nbasis); % prior mean and variance for basis coefficients of template function
prior_params.alpha = 4; prior_params.beta = 0.01; % hyperparameters for \sigma^2
stepvar_chain_ini = 0.05*cov(((q_basis'*q_basis)\q_basis'*q_mat)'); % initial proposal variance of 'c'
checknum_stop_fix = 1000; % the number of proposal variance adjustment
div_iter_by = 1000; % proposal variance adjustment every 'div_iter_by' number of MCMC iterations
target_rate = [0.2 0.55]; % target MH acceptance rate

qs = q_mat(:,1:N_func); % given SRVFs
total_mcmc = 10000; % total number of MCMC iterations

% run MCMC
[c_samples, d_samples, sigma_samples, log_post, kappa_prop, ...
    MCMC_time] = mcmc_regist(total_mcmc, c_ini, ds_ini, ...
    sigma_squared_ini, t, qs, q_basis, prior_params, kappa_prop, ...
    stepvar_chain_ini, target_rate, checknum_stop_fix,div_iter_by);

%% Sequential Bayesian Updating
% Assign equal weights to posterior samples obtained using MCMC
% Update weighted posterior samples through SMC when a new function arrives

J_particles = 10000; % number of posterior particles

% burn-in MCMC posterior samples and assign equal weights
c = c_samples((end-J_particles+1):end,:);
ds = d_samples((end-J_particles+1):end,:,:);
sigma_squared = sigma_samples((end-J_particles+1):end);
weights = 1/J_particles*ones(1,J_particles);

% proposal concentration parameter for phase components
kappa_prop = [kappa_prop(1:N_func);mean(kappa_prop(1:N_func))*ones(N-N_func,1)];

c_smc = cell(N,1);
d_smc = cell(N,1);
sigma_squared_smc = zeros(J_particles,N);
w_smc = zeros(J_particles,N);
ESS_smc = zeros(2,N);
logpost_smc = zeros(J_particles,N);
SMC_timerecord = zeros(N,1);

K_mh_c = 30;
K_mh_d = 30;

for N_func = (N_func+1):N

    qs = q_mat(:,1:N_func); % given SRVFs

    % run SMC algorithm
    [c_update, d_update, sigma_squared_update, w_update,...
        ESS, logpost_update, SMC_time] = smc_regist(c, ds, ...
        sigma_squared, weights, t, qs, q_basis, ...
        prior_params, kappa_prop, K_mh_c, K_mh_d);

    c_smc{N_func} = c_update;
    d_smc{N_func} = d_update;
    sigma_squared_smc(:,N_func) = sigma_squared_update;
    w_smc(:,N_func) = w_update;
    ESS_smc(:,N_func) = ESS;
    logpost_smc(:,N_func) = logpost_update;
    SMC_timerecord(N_func) = SMC_time;
    
    c = c_update;
    ds = d_update;
    sigma_squared = sigma_squared_update;
    weights = w_update;
end

