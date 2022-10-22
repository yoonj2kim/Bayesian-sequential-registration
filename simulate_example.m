% setup paths
addpath(genpath('DP'))
addpath(genpath('bspline_tools'))
addpath(genpath('basis'))

%% Data generation

M = 100; % the number of evalutated x values
t = linspace(0,1,M); % the x values where functions are evaluated at

N = 100; % the number of functions

%%%%% Template function %%%%%

nbasis = 8; % the number of basis for the template function
q_basis = basis_bspline(t, nbasis, 4).matrix; % create a B-spline design matrix
c_mu = [1.0120;2.4743;2.0992;-4.1448;4.1974;-2.3106;-2.0217;-1.8065]; % basis coefficients for a two-peaked template function
q_mu = q_basis*c_mu; % template SRVF

%%%%% Add randomly drawn phase components and noises %%%%%

q_mat = zeros(M,N); % generated SRVFs
f_mat = zeros(M,N); % generated functions

kappa = 50;  % concentration parameter for phase component generation
M_gamma = 5; tG = linspace(0,1,M_gamma)'; % phase component partition size
gamma = zeros(N,M); % phase functions evaluated at 't'
gamma_M = zeros(N,M_gamma); % phase functions evaluated at 'tG'
mat = (diag(ones(M_gamma-1,1),1)+diag(-ones(M_gamma,1))); % matrix to derive phases using Dirichlet process
mat = mat(1:(end-1),:);  
sigma_squared_gt = 0.03; % variance of error process

for i = 1:N
    % generate phase increments at 'tG' using Dirichlet process
    increments_star = gamrnd(mat*tG*kappa,1);
    increments_star = increments_star./sum(increments_star);
    gamma_M(i,:) = [0;cumsum(increments_star)];
    gamma_M(i,end) = 1;
    % evaluate phase components at 't'
	d_temp = interp1(tG,gamma_M(i,:),t,'linear'); 
	d_temp(end) = 1;
	gamma(i,:) = d_temp;
    % generated SRVFs: N((q_\mu, \gamma_i^{-1})([t]), \sigma^2 I_M)
    q_mat(:,i) = warp_q_gamma(q_mu,invertGamma(gamma(i,:)),t) + ...
        mvnrnd(zeros(M,1),sigma_squared_gt);
    % generated functions at the original function space
    f_mat(:,i) = srvf_to_f(q_mat(:,i),t,0);
end

%% Initialization
% Derive posterior samples given functions 1 through
% N_func using MCMC.

N_func = 4;

% initial samples
c_ini = zeros(1,nbasis);
ds_ini = repmat(tG',[N_func 1]);
sigma_squared_ini = 0.1;

% set hyperparameters and tuning parameters
kappa_prop = kappa*ones(N_func,1); % proposal concentration parameter for phase components
prior_params.kappa = 5; % prior concentration parameter for phase components
prior_params.c_mean = zeros(nbasis,1);
prior_params.c_var = 20*eye(nbasis); % prior mean and variance for basis coefficients of template function
prior_params.alpha = 4; prior_params.beta = 0.01; % hyperparameters for \sigma^2
stepvar_chain_ini = 0.05*cov(((q_basis'*q_basis)\q_basis'*q_mat)'); % initial proposal variance of 'c'
checknum_stop_fix = 1000; % the number of proposal variance adjustment
div_iter_by = 1000; % proposal variance adjustment every 'div_iter_by' number of MCMC iterations
target_rate = [0.2 0.55]; % target MH acceptance rate

qs = q_mat(:,1:N_func); % given SRVFs
total_mcmc = 500000; % total number of MCMC iterations

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

%% Visualization of Target Posteriors

% target posterior given N_func functions
N_func = 6;
c_samples = c_smc{N_func};
d_samples = d_smc{N_func};
sigma_squared_samples = sigma_squared_smc(:,N_func);
weights = w_smc(:,N_func);

weight_scale = 20; % scale up weights for visualization of posterior uncertainty
N_func_plot = 5; % # of functions to plot

N_samples_plot = 1000; % # of posterior samples plotted to show posterior uncertainty
J_particles = length(weights);
index_plot = round(linspace(1,J_particles,N_samples_plot));
lw = 2; % linewidth

% axis limits and figure settings
minf = NaN;
maxf = NaN;
for i = 1:N_func
    minf = min([minf min(f_mat(:,i))]);
    maxf = max([maxf max(f_mat(:,i))]);
end
ylim_plot = [minf-(maxf-minf)/6 maxf+(maxf-minf)/6];
ytick_plot = [0 0.5];
xtick_plot = [0 1];
font_size = 12;
plot_color = 'black';

% plot given functions
for i = 1:N_func_plot
    figure
    plot(t,f_mat(:,i),'Color',plot_color,'linewidth',lw)
    ylim(ylim_plot)
    title(sprintf('$f_{%d}$',i),'Interpreter','latex','Fontsize',font_size)
    set(gca,'xtick',xtick_plot,'ytick',ytick_plot,'FontSize',font_size)
end

% marginal posterior uncertainty of phase components
for i = 1:N_func_plot
    figure
    for j = 1:N_samples_plot
        p = plot(tG,d_samples(index_plot(j),:,i),...
                'Color',plot_color,'linewidth',lw);    
        p.Color(4) = weight_scale*weights(index_plot(j));        
            hold on
    end    
    plot(tG,gamma_M(i,:),'r','linewidth',lw/2)
    axis square
    title(strjoin({'$\gamma_{',num2str(i),'}$'}),'Interpreter','latex','Fontsize',font_size)
    set(gca,'xtick',xtick_plot,'ytick',xtick_plot,'FontSize',font_size)
end

% given functions after applying the corresponding phase posterior samples
for i = 1:N_func_plot
    figure
    hold off
    for j = 1:N_samples_plot
        d_temp = d_samples(index_plot(j),:,i);
        d_temp = interp1(tG,d_temp,t','linear');
        p = plot(t,warp_f_gamma(f_mat(:,i),d_temp,t),...
            'Color',plot_color,'linewidth',lw);   
        p.Color(4) = weight_scale*weights(index_plot(j));        
            hold on
    end    
    ylim(ylim_plot)
    title(strjoin({'$f_{',num2str(i),'}\circ\gamma_{',num2str(i),'}$'}),'Interpreter','latex','Fontsize',font_size)
    set(gca,'xtick',xtick_plot,'FontSize',font_size)
    set(gca,'ytick',ytick_plot,'FontSize',font_size)
end

% summarize posterior of the template function
[f_mean,q_mean,f_var_ptw,q_var_ptw,eFR_var,f_samples] = template_posterior(c_samples,q_basis,weights,mean(f_mat(1,1:N_func)));

figure
hold on
for j = 1:N_samples_plot
    p = plot(t,f_samples(:,index_plot(j)),'k','linewidth',lw);
    p.Color(4) = weight_scale*weights(index_plot(j));         
end    
plot(t,srvf_to_f(q_mu,t,0),'r','linewidth',lw/2);
ylim(ylim_plot)
title('Marginal posterior of template function','Fontsize',font_size)
set(gca,'xtick',xtick_plot,'ytick',ytick_plot,'FontSize',font_size)
 
% marginal posterior uncertainty of sigma squared
figure
hold on
sigma_squared_temp = linspace(0,0.05,100);
marg_post_pdf = exp(log_invgamma_pdf(sigma_squared_temp,prior_params.alpha,prior_params.beta));
plot(sigma_squared_temp,marg_post_pdf/max(marg_post_pdf),'linewidth',lw,'Color','b')
[marg_post_pdf,sigma_squared_temp] = ksdensity(sigma_squared_samples, ...
    'Support',[sigma_squared_temp(1) sigma_squared_temp(end)],'Weights',weights);
plot(sigma_squared_temp,marg_post_pdf/max(marg_post_pdf),'linewidth',lw,'Color','k')
xline(sigma_squared_gt,'Color','r','LineWidth',1.5)
xlim([0 0.05])
set(gca,'FontSize',font_size)



