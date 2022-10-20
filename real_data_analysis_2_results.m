%% load posterior particles given El Nino years

load('data/real_data_analysis_2_elnino_data.mat')

c_samples = c_smc{N_func};
d_samples = d_smc{N_func};
weights = w_smc(:,N_func);

[f_mean_el,q_mean_el,f_var_ptw_el,q_var_ptw_el,eFR_var_el, ...
    f_samples_el] = template_posterior(c_samples,q_basis,weights,mean(f_mat(1,:)));

xconf = [t t(end:-1:1)];
yconf_el = [f_mean_el+2*sqrt(f_var_ptw_el') f_mean_el(end:-1:1)-2*sqrt(f_var_ptw_el(end:-1:1)')];
ft_el = f_mat;

[d_mean_el,~] = phase_posterior(d_samples, weights);

%% load posterior particles given La Nina years

load('real_data_analysis_2_lanina_data.mat')

c_samples = c_smc{N_func};
d_samples = d_smc{N_func};
weights = w_smc(:,N_func);

[f_mean_la,q_mean_la,f_var_ptw_la,q_var_ptw_la,eFR_var_la, ...
    f_samples_la] = template_posterior(c_samples,q_basis,weights,mean(f_mat(1,:)));

yconf_la = [f_mean_la+2*sqrt(f_var_ptw_la') f_mean_la(end:-1:1)-2*sqrt(f_var_ptw_la(end:-1:1)')];
ft_la = f_mat;

[d_mean_la,~] = phase_posterior(d_samples, weights);

%% load posterior particles given el nino years

load('real_data_analysis_2_neutral_data.mat')

[f_mean_neu,q_mean_neu,f_var_ptw_neu,q_var_ptw_neu,eFR_var_neu, ...
    f_samples_neu] = template_posterior(c_samples,q_basis,weights,mean(f_mat(1,:)));

yconf_neu = [f_mean_neu+2*sqrt(f_var_ptw_neu') f_mean_neu(end:-1:1)-2*sqrt(f_var_ptw_neu(end:-1:1)')];
ft_neu = f_mat;

[d_mean_neu,~] = phase_posterior(d_samples, weights);

%% posterior uncertainties of the template function 

figure

p = fill(xconf,yconf_el,'red');
p.FaceColor = [1 0 0];      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.09;
hold on
plot(t,f_mean_el,'Color',[1 0 0],'LineWidth',1.5)

p = fill(xconf,yconf_la,'yellow');
p.FaceColor = [0.9290 0.6940 0.1250];      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.2;
hold on
plot(t,f_mean_la,'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)

p = fill(xconf,yconf_neu,'blue');
p.FaceColor = [0 0.4470 0.7410];      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.15;
hold on
plot(t,f_mean_neu,'Color',[0 0.4470 0.7410],'LineWidth',1.5)

%% summary of registration of El Nino functions using posterior mean of phase components

figure
for i = 1:size(ft_el,2)
    plot(t,ft_el(:,i),'linewidth',1.5)
    hold on
end

figure
for i = 1:size(ft_el,2)
    plot(tG,d_mean_el(:,i),'linewidth',1.5)
    hold on
end

figure
for i = 1:size(ft_el,2)
    plot(t,warp_f_gamma(ft_el(:,i),interp1(tG,d_mean_el(:,i),t),t),'linewidth',1.5)
    hold on
end

%% summary of registration of La Nina functions using posterior mean of phase components

figure
for i = 1:size(ft_la,2)
    plot(t,ft_la(:,i),'linewidth',1.5)
    hold on
end

figure
for i = 1:size(ft_la,2)
    plot(tG,d_mean_la(:,i),'linewidth',1.5)
    hold on
end

figure
for i = 1:size(ft_la,2)
    plot(t,warp_f_gamma(ft_la(:,i),interp1(tG,d_mean_la(:,i),t),t),'linewidth',1.5)
    hold on
end


%% summary of registration of neutral years using posterior mean of phase components

figure
for i = 1:size(ft_neu,2)
    plot(t,ft_neu(:,i),'linewidth',1.5)
    hold on
end

figure
for i = 1:size(ft_neu,2)
    plot(tG,d_mean_neu(:,i),'linewidth',1.5)
    hold on
end

figure
for i = 1:size(ft_neu,2)
    plot(t,warp_f_gamma(ft_neu(:,i),interp1(tG,d_mean_neu(:,i),t),t),'linewidth',1.5)
    hold on
end

