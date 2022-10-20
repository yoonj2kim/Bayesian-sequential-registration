
load('data/simulated_example_2_data.mat')

c_samples = c_smc{N_func};
d_samples = d_smc{N_func};
sigma_squared_samples = sigma_squared_smc(:,N_func);
weights = w_smc(:,N_func);

weight_scale = 20; % scale up weights for visualization of posterior uncertainty
N_func_plot = N; % # of functions to plot

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

[d_mean,d_var_ptw,post_var] = phase_posterior(d_samples,weights);

% multimodality of 7th phase component
d_temp = squeeze(d_samples(:,:,7));

mode_1_idx = [];
mode_2_idx = [];
for j = 1:J_particles
    if d_temp(j,3) > tG(3)
        mode_1_idx = [mode_1_idx j];
    else
        mode_2_idx = [mode_2_idx j];
    end
end

d7_mean_mode_1 = d_temp(mode_1_idx,:)'*weights(mode_1_idx)/sum(weights(mode_1_idx));
d7_mean_mode_2 = d_temp(mode_2_idx,:)'*weights(mode_2_idx)/sum(weights(mode_2_idx));

% marginal posterior uncertainty of phase components
for i = 1:6
    figure
    plot(tG, d_mean(:,i))
    axis square
    title(strjoin({'$\gamma_{',num2str(i),'}$'}),'Interpreter','latex','Fontsize',font_size)
    set(gca,'xtick',xtick_plot,'ytick',xtick_plot,'FontSize',font_size)
end

figure
i = 7;
plot(tG,d7_mean_mode_1,'Color','red','linewidth',1);
hold on
plot(tG,d7_mean_mode_2,'Color','blue','linewidth',1);
axis square
ylim([0 1])
title(strjoin({'$\gamma_{',num2str(i),'}$'}),'Interpreter','latex','Fontsize',font_size)
set(gca,'xtick',xtick_plot,'ytick',xtick_plot,'FontSize',font_size)

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
    set(gca,'xtick',xtick_plot,'ytick',ytick_plot,'FontSize',font_size)
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
 
