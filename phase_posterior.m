function [d_mean,d_var_ptw,post_var] = phase_posterior(d_samples,weights)

    [N_samples,nG,N] = size(d_samples);
    tG = linspace(0,1,nG);
    
    d_mean = zeros(nG,N);
    for i = 1:N
        d_mean(:,i) = squeeze(d_samples(:,:,i))'*weights;
    end
    
    d_var_ptw = zeros(nG,N);
    for i = 1:N
        for j = 2:(nG-1)
            pw_var_temp = 0;
            for s = N_samples
                pw_var_temp = pw_var_temp + weights(s)*(d_samples(s,j,i)-d_mean(j,i))^2;
            end
            d_var_ptw(j,i) = pw_var_temp;
        end
    end
    
    post_var = zeros(N,1);
    for i = 1:N
        var_temp = 0;
        for s = N_samples
            var_temp = var_temp + weights(s)*trapz(tG,(squeeze(d_samples(s,:,i))'-d_mean(:,i)).^2);
        end
        post_var(i) = var_temp;
    end
end