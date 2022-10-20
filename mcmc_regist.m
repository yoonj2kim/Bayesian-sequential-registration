function [c_samples, d_samples, sigma_squared_samples, log_post, kappa_prop, ...
    MCMC_time] = mcmc_regist(iterations, c_ini, ds_ini, ...
    sigma_squared_ini, t, qs, q_basis, prior_params, kappa_prop, ...
    stepvar_chain_ini, target_rate, checknum_stop_fix, div_iter_by)
    
    %% set up
    % check dimensions of given matrices for further calculations
    if size(t,2) ~= 1 && size(t,1) == 1
        t = t'; 
    end
    
    if size(qs,1) ~= length(t)
        qs = qs';
    end
    
    if size(qs,2) ~= size(ds_ini,2)
        ds_ini = ds_ini';
    end

    if size(c_ini,1) ~= 1
        c_ini = c_ini';
    end
    
    N_func = size(qs,2);
    nbasis = length(c_ini); % the number of basis for the template function
    M_gamma = size(ds_ini,1); tG = linspace(0,1,M_gamma)'; % phase component partition size
    mat = (diag(ones(M_gamma-1,1),1)+diag(-ones(M_gamma,1))); % matrix to derive phases using Dirichlet process
    mat = mat(1:(end-1),:);  
    
    % set hyperparameters and tuning parameters
    kappa_prior = prior_params.kappa;
    c_mean = prior_params.c_mean;
    c_var = prior_params.c_var;
    alpha_sigma = prior_params.alpha;
    beta_sigma = prior_params.beta;

    stepvar_chain = stepvar_chain_ini;
    temp_stepvar_chain = stepvar_chain;
    naccept_c = zeros(checknum_stop_fix,1);
    naccept_d = zeros(checknum_stop_fix,N_func);
    rate = zeros(checknum_stop_fix,1);

    log_post = zeros(iterations,1); % unnormalized log posterior
    d_samples = zeros(iterations,M_gamma,N_func);
    c_samples = zeros(iterations,nbasis);
    sigma_squared_samples = zeros(iterations,1);
 
    % initialization

    sigma_squared_curr = sigma_squared_ini;
    d_curr = ds_ini;
    c_curr = c_ini;    

    q_mu_curr = c_curr*q_basis';
    gam_curr = interp1(tG,d_curr,t,'linear');
    SSE_curr = SSE_q(gam_curr,t,qs,q_mu_curr);
    logl_cur = f_logl(gam_curr,t, qs, q_mu_curr, sigma_squared_curr, SSE_curr); % log likelihood
    prior_c_cur = logmvnpdf(c_curr,c_mean',c_var); % log prior density of template basis coefficients
    
    % log prior density of phase components
    l_pi_cur = zeros(N_func,1);   
    for ii = 1:N_func
        increments_cur = mat*d_curr(:,ii);
        l_pi_cur(ii) = sum(log(increments_cur + eps).*(kappa_prior*mat*tG));
    end
    
    k=1;   

    %% Run MCMC
    tic

    for i = 1:iterations
        %%%%% template function MH update %%%%%
        stepvar_chain = 0.5*(stepvar_chain + stepvar_chain');
        [~,p]=chol(stepvar_chain);
        if p ~=0
            stepvar_chain = diag(diag(stepvar_chain));
        end
        
        c_prop = mvnrnd(c_curr, stepvar_chain);
        q_star_prop = c_prop*q_basis'; 
        SSE_prop = SSE_q(gam_curr,t, qs, q_star_prop);
        logl_prop = f_logl(gam_curr,t, qs, q_star_prop, sigma_squared_curr, SSE_prop);
        prior_c_prop = logmvnpdf(c_prop,c_mean',c_var);
        
        log_accept_prob = logl_prop+prior_c_prop-logl_cur-prior_c_cur;

        if (log(rand) < log_accept_prob)
            SSE_curr = SSE_prop;
            q_mu_curr = q_star_prop;
            naccept_c(k,:) = naccept_c(k,:) + 1;
        end
        
        %%%%% gamma components MH update %%%%%
        % tune proposal concentration parameter
        for ii = 1:N_func  
        if(mod(i,div_iter_by) == 0)
            if((naccept_d(k,ii)/div_iter_by)>target_rate(2) ||(naccept_d(k,ii)/div_iter_by)<target_rate(1)) 
                kappa_prop(ii) = max(10,kappa_prop(ii)*.39/(naccept_d(k,ii)/div_iter_by + 1e-10));
            end
        end
        
        % propose candidate gamma component 
        increments_star = gamrnd(mat*tG*kappa_prop(ii),1);
        increments_star = increments_star./sum(increments_star);
        gamma_star = [0;cumsum(increments_star)];
        d_prop = interp1(tG,d_curr(:,ii),gamma_star,'linear');
        d_prop(end) = 1;
        increments_can = mat*d_prop;
        
        % evaluate log prior density
        l_pi_can = sum(log(increments_can + eps).*(kappa_prior*mat*tG));
        
        % evaluate likelihoods
        gam_prop = interp1(tG,d_prop,t,'linear');
        gam_prop(end) = 1;
        SSE_prop = SSE_q(gam_prop,t,qs(:,ii),q_mu_curr);
        star_inverse = interp1(gamma_star,tG,tG,'linear');
        
        % evaluate proposal densities
        star_inverse(end) = 1;
        increments_star_inverse = mat*star_inverse;
        increments_star = mat*gamma_star;
        l_g_cur_given_can = sum(log(increments_star_inverse + eps).*(kappa_prop(ii)*mat*tG));
        l_g_can_given_cur = sum(log(increments_star + eps).*(kappa_prop(ii)*mat*tG));
        
        % evaluate Jacobian 
        grad_gamma_cur_inv = gradient(interp1(tG,tG,d_curr(:,ii)),tG);
        grad_gamma_can_inv = gradient(interp1(tG,tG,d_prop),tG);
        jacobian = sum(log(interp1(tG,grad_gamma_can_inv,d_curr(2:(end-1),ii))))-sum(log(interp1(tG,grad_gamma_cur_inv,d_prop(2:(end-1)))));

        % evaluate Acceptance Probability
        log_accept_prob = (SSE_curr(ii)-SSE_prop)/(2*sigma_squared_curr) + (l_pi_can+l_g_cur_given_can)-(l_pi_cur(ii)+l_g_can_given_cur) + jacobian;
        
        % evaluate Acceptance
        if(log(rand) < log_accept_prob)
            d_curr(:,ii) = d_prop;
            naccept_d(k,ii) = naccept_d(k,ii) + 1;
            SSE_curr(ii) = SSE_prop;
            l_pi_cur(ii) = l_pi_can;
        end    
        end
        
        %%%% Centering Step %%%%%
        gamma_mean = SqrtMeanInverse(d_curr'); 
        for ii = 1:N_func
            d_curr(:,ii) = interp1(tG, d_curr(:,ii)', (tG(end)-tG(1)).*gamma_mean + tG(1));
            increments_cur = mat*d_curr(:,ii);
            l_pi_cur(ii) = sum(log(increments_cur + eps).*(kappa_prior*mat*tG));
        end
        
        q_mu_curr = warp_q_gamma(q_mu_curr,interp1(tG,gamma_mean,t),t);
        c_curr = q_mu_curr*q_basis/(q_basis'*q_basis);
        q_mu_curr = c_curr*q_basis';     
        
        prior_c_cur = logmvnpdf(c_curr,c_mean',c_var);
        gam_curr = interp1(tG,d_curr,t,'linear');
        SSE_curr = SSE_q(gam_curr,t,qs,q_mu_curr);
        logl_cur = f_logl(gam_curr,t, qs, q_mu_curr, sigma_squared_curr, SSE_curr);

        %%%%% Gibbs update for sigma_squared %%%%%
        alpha_conditional = alpha_sigma + length(t)*N_func/2;
        beta_conditional = beta_sigma + sum(SSE_curr)/2;
        invS = gamrnd(alpha_conditional,1/beta_conditional);
        sigma_squared_curr = 1/invS;
        
        %%%%% Record MCMC samples and unnormalized log posterior %%%%%
        log_post(i) = prior_c_cur + sum(logl_cur) + sum(l_pi_cur) + log_invgamma_pdf(sigma_squared_curr,alpha_sigma,beta_sigma);    
        d_samples(i,:,:) = d_curr;
        c_samples(i,:) = c_curr;
        sigma_squared_samples(i) = sigma_squared_curr;

    %%%%% Tune proposal variance for template basis coefficients %%%%%
    if rem(i,div_iter_by) == 0
        k = k + 1;
        checknum_ind = i/div_iter_by;
        rate(checknum_ind,:) = ( 1 + naccept_c(checknum_ind,:))./(div_iter_by + 2);
            for aa = 1:nbasis
                if (checknum_ind <= checknum_stop_fix) && (rate(checknum_ind,:) < target_rate(1) || rate(checknum_ind,:) > target_rate(2))  % Note that I stop making adjustments by kjk_stop_fix iterations
                    stepvar_chain(aa,aa) = stepvar_chain(aa,aa)*rate(checknum_ind,:)/(target_rate(1) + range(target_rate)/2);
                end
            end
        if checknum_ind <= checknum_stop_fix
        try
            span = (i-div_iter_by+1):i;
            chaincorr = corr(c_samples(span,:));
            chaincorr(isnan(chaincorr)) = 0;
            new_step_cor = zeros(nbasis,nbasis);
                  for rr = 1:nbasis
                        for cc = 1:nbasis
                            if (rr==cc)
                                new_step_cor(rr,cc) = 1;
                            elseif (abs(chaincorr(rr,cc))<0.4)
                                new_step_cor(rr,cc)=0;
                            else
                                new_step_cor(rr,cc) = sign(chaincorr(rr,cc))*(abs(chaincorr(rr,cc))-0.2);
                            end
                        end
                  end
                  for rr = 1:nbasis
                        for cc = 1:nbasis
                            stepvar_chain(rr,cc) = new_step_cor(rr,cc)*sqrt(stepvar_chain(rr,rr)*stepvar_chain(cc,cc));
                        end
                  end
                  def_by = 10;
                  pvars = diag(cov(c_samples(span,:)))/def_by;
                    for rr = 1:nbasis
                        for cc = 1:nbasis
                            temp_stepvar_chain(rr,cc) = new_step_cor(rr,cc)*sqrt(pvars(rr)*pvars(cc));
                        end
                    end
                    tempvar = diag(temp_stepvar_chain);
                    for dd = 1:nbasis
                        if tempvar(dd)<=0
                            temp_stepvar_chain(dd,dd) = stepvar_chain(dd,dd);
                        end
                    end
                    [~,posdef] = chol(temp_stepvar_chain);
                        if posdef ~= 0 || ~isnan(sum(sum(temp_stepvar_chain)))
                            stepvar_chain = temp_stepvar_chain;
                        end
        catch
            disp('Error in covariance update')
        end
        end
    end
    end
    
    MCMC_time = toc;
end

function [logp] = logmvnpdf(x,mu,Sigma)
    % log pdf of MVN(mu,Sigma)
    
    [~,D] = size(x);
    const = -0.5 * D * log(2*pi);
    
    xc = bsxfun(@minus,x,mu);
    [~,p]=chol(Sigma);
    if p ==0
        term1 = -0.5 * sum((xc / Sigma) .* xc, 2);
        term2 = const - 0.5 * logdet(Sigma);
        logp = term1' + term2;
    else
        logp = -inf;
    end

end

function y = logdet(A)

    U = chol(A);
    y = 2*sum(log(diag(U)));

end