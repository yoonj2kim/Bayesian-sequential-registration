function [f_mean,q_mean,f_var_ptw,q_var_ptw,eFR_var,f_samples] = template_posterior(c_samples,...
    q_basis,weights,f0)
    M = size(q_basis,1);
    N_samples = size(c_samples,1);
    
    t = linspace(0,1,M);
    q_samples = q_basis*c_samples';
    q_mean = sum(q_samples'.*weights);
    f_samples = zeros(M,N_samples);
    
    for j = 1:N_samples
        f_samples(:,j) = srvf_to_f(q_basis*c_samples(j,:)',t,f0);
    end
    
    f_mean = sum(f_samples'.*weights);
    f1_min = min(f_samples(end,:));
    f1_max = max(f_samples(end,:));
    f1_mean = f_mean(end);
    f1_range = f1_max-f1_min;

    b = f1_mean+f1_range/4;
    a = f1_mean-f1_range/4;

    for j = 1:N_samples
        f_samples(:,j) = f_samples(:,j)-f_samples(end,j)+(a+(f_samples(end,j)-f1_min)*(b-a)/(f1_max-f1_min));
    end  
    
    f_mean = sum(f_samples'.*weights);
    f_var_ptw = zeros(M,1);
    for i = 1:M
        pw_var_temp = 0;
        for j = 1:N_samples
            pw_var_temp = pw_var_temp + weights(j)*(f_samples(i,j)-f_mean(i))^2;
        end
        f_var_ptw(i) = pw_var_temp;
    end

    q_var_ptw = zeros(M,1);
    for i = 1:M
        pw_var_temp = 0;
        for j = 1:N_samples
            pw_var_temp = pw_var_temp + weights(j)*(q_samples(i,j)-q_mean(i))^2;
        end
        q_var_ptw(i) = pw_var_temp;
    end
    
    var_temp = 0;
    for j = 1:N_samples
        var_temp = var_temp + weights(j)*trapz(t,(q_samples(:,j)-q_mean').^2);
    end

    eFR_var = var_temp;
end