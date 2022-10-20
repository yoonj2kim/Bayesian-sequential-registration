
load('simulated_example_1_data.mat')


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
plot(t,srvf_to_f(q_mu,t,0),'r','linewidth',lw/2);
ylim(ylim_plot)
title('Marginal posterior of template function','Fontsize',font_size)
set(gca,'xtick',xtick_plot,'ytick',ytick_plot,'FontSize',font_size)
 
% figure with pointwise standard deviation bands around the mean to enhance the display
max_temp_var = max(sqrt(f_var_ptw));
min_temp_var = min(sqrt(f_var_ptw));
scale_up_sd = 20;

figure
patch([t NaN],[f_mean NaN],[sqrt(f_var_ptw)' NaN],...
    'EdgeColor','interp',...
    'MarkerFaceColor','flat',...
    'LineWidth',lw);
colormap turbo;
%colorbar;
clim([min_temp_var max_temp_var])
hold on
plot(t,f_mean+scale_up_sd*sqrt(f_var_ptw'),'Color',[0.6 0.6 0.6])
plot(t,f_mean-scale_up_sd*sqrt(f_var_ptw'),'Color',[0.6 0.6 0.6])
ylim(ylim_plot)
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


% ESS
figure
plot(31:100,ESS_smc(2,31:100),'.','MarkerSize', 10)
line([30 100],[5000 5000],'LineStyle','--','Color','black')
ylim([0 10000])
set(gca,'xtick',[30 60 90],'FontSize',12)
set(gca,'ytick',[0 5000 10000],'FontSize',12)

