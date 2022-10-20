
load('real_data_analysis_3_data.mat')

N_func = 50;

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
ytick_plot = [0 1];
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
    axis square
    title(strjoin({'$\gamma_{',num2str(i),'}$'}),'Interpreter','latex','Fontsize',font_size)
    set(gca,'xtick',xtick_plot,'ytick',xtick_plot,'FontSize',font_size)
end

% summarize posterior of the template function
[f_mean,q_mean,f_var_ptw,q_var_ptw,eFR_var,f_samples] = template_posterior(c_samples,q_basis,weights,mean(f_mat(1,1:N_func)));

figure
hold on
for j = 1:N_samples_plot
    p = plot(t,f_samples(:,index_plot(j)),'k','linewidth',lw);
    p.Color(4) = weight_scale*weights(index_plot(j));         
end 
ylim(ylim_plot)
title('Marginal posterior of template function','Fontsize',font_size)
set(gca,'xtick',xtick_plot,'ytick',ytick_plot,'FontSize',font_size)

[d_mean,d_var_ptw,post_var] = phase_posterior(d_samples,weights);

% marginal posterior means of phase components
figure
hold on
for i = 1:N_func
    plot(tG,d_mean(:,i),'linewidth',lw)
end

% given functions registered using marginal posterior means of phase components
figure
hold on
for i = 1:N_func
    gamma_temp = interp1(tG,d_mean(:,i),t','linear');
    plot(t,warp_f_gamma(f_mat(:,i),gamma_temp,t),'linewidth',lw)
end

 