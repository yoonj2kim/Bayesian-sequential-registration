function [c_update, d_update, sigma_squared_update, w_update, ESS, logpost_update, ...
    SMC_time] = smc_regist(c, ds, sigma_squared, weights, t, qs, q_basis, prior_params, ...
    kappa_prop, K_mh_c, K_mh_d)
    
    %% set up
    % check dimensions of given matrices for further calculations
    if size(t,2) ~= 1 && size(t,1) == 1
        t = t'; 
    end
    
    if size(qs,1) ~= length(t)
        qs = qs';
    end

    if size(c,1) ~= length(weights)
        c = c';
    end
    
    if size(weights,2) == 1
        weights = weights';
    end

    N_func = size(qs,2); % # of functions
    nbasis = size(c,2); % # of basis functions for the template function
    J_particles = length(weights); % # of particles
    M_gamma = size(ds,2); tG = linspace(0,1,M_gamma)'; % phase component partition size
    mat = (diag(ones(M_gamma-1,1),1)+diag(-ones(M_gamma,1))); % matrix to derive phases using Dirichlet process
    mat = mat(1:(end-1),:);  

    % set hyperparameters and tuning parameters
    kappa_prior = prior_params.kappa;
    c_mean = prior_params.c_mean;
    c_var = prior_params.c_var;
    alpha_sigma = prior_params.alpha;
    beta_sigma = prior_params.beta;

    logpost_update = zeros(J_particles,1);
    q_new = qs(:,end);

    %% Run SMC
    tic

    % Initialize the new phase component using optimization
    
    w_samples = zeros(J_particles,1);
    q_mu_samples = c*q_basis';
    d_samples = ds;
    d_samples(:,:,N_func) = zeros(J_particles,M_gamma);
    c_samples = c;

    parfor i = 1:J_particles
        q_mu_curr = q_mu_samples(i,:)';   
        gam = optimum_reparam(q_mu_curr,q_new,t); % find optimal phase function
        gam_ini = lsq_lut_piecewise(t,gam,tG); % replace with a piecewise linear least squares fit

        % make sure phase components are non-decreasing
        increments_cur = mat*gam_ini;
        nondecreasing_check = sum(increments_cur<0);
        gam_ini(1) = 0;
        gam_ini(end) = 1;
        
        while nondecreasing_check ~= 0
            ind_neg = find(increments_cur<0);
            for j = 1:length(ind_neg)
                if ind_neg(j) < M_gamma-1
                    gam_ini(ind_neg(j)+1) = 0.5*(gam_ini(ind_neg(j))+gam_ini(ind_neg(j)+2));
                else
                    gam_ini(ind_neg(j)) = 0.5*(gam_ini(ind_neg(j)-1)+gam_ini(ind_neg(j)+1));
                end
            end
            
            gam_ini(1) = 0;
            gam_ini(end) = 1;
            increments_cur = mat*gam_ini;
            nondecreasing_check = sum(increments_cur<0);
        end

        gam_ini(1) = 0;
        gam_ini(end) = 1;
        gam_ini_interp = interp1(tG,gam_ini,t,'linear');

        increments_cur = mat*gam_ini;
        l_pi_cur_new = sum(log(increments_cur + eps).*(kappa_prior*mat*tG));
        SSE_ini = SSE_q(gam_ini_interp, t, q_new, q_mu_curr);
        logl_ini = f_logl(gam_ini_interp, t, q_new, q_mu_curr, sigma_squared(i), SSE_ini);
        w_samples(i) = weights(i)*exp(l_pi_cur_new + logl_ini);
        d_temp = [squeeze(ds(i,:,:)) gam_ini];
        d_samples(i,:,:) = d_temp;
    end    

    % check ESS and resample if necessary 
    if sum(w_samples) == 0
        w_samples = ones(J_particles,1);
    end
    w_samples = w_samples/sum(w_samples);
    ESS = 1/sum(w_samples.^2);
    if ESS < J_particles/2
        resam_num = mnrnd(J_particles,w_samples);
        resam_idx = cumsum(resam_num);
        c_resam = zeros(size(c_samples));
        d_resam = zeros(size(d_samples));
        sigma_resam = zeros(length(sigma_squared),1);
        for i = 1:J_particles
            if resam_num(i) > 0
                if i == 1
                    start_idx = 1;
                else
                    start_idx = resam_idx(i-1)+1;
                end
                c_resam(start_idx:(start_idx+resam_num(i)-1),:) = repmat(c_samples(i,:), [resam_num(i) 1]);
                d_resam(start_idx:(start_idx+resam_num(i)-1),:,:) = repmat(d_samples(i,:,:), [resam_num(i) 1]);
                sigma_resam(start_idx:(start_idx+resam_num(i)-1)) = repmat(sigma_squared(i), [resam_num(i) 1]);
            end
        end
        d_samples = d_resam;
        sigma_squared = sigma_resam;
        c_samples = c_resam;
        w_samples = repmat(1/J_particles,[J_particles 1]);
    end

    % tune proposal variance for the template function
    stepmean_chain = sum(c_samples.*repmat(w_samples,[1 nbasis]),1);
    stepvar_chain = zeros(nbasis);
    for i = 1:J_particles
        stepvar_chain = stepvar_chain + w_samples(i)*...
            (c_samples(i,:)-stepmean_chain)'*(c_samples(i,:)-stepmean_chain);
    end
    stepvar_chain = 0.5*(stepvar_chain+stepvar_chain');
    [~,p] = chol(stepvar_chain);
    if p ~= 0
        stepvar_chain = diag(diag(stepvar_chain));
    end

    %%%%% Perturb augmented particles %%%%%
    sigma_samples = sigma_squared;
    parfor i = 1:J_particles
        % set the initial sample as the ith sample in the previous sequence
        sigma_cur = sigma_squared(i);
        d_curr = squeeze(d_samples(i,:,:));
        c_curr = c_samples(i,:);
        q_mu_curr = c_curr*q_basis';      
        prior_c_curr = logmvnpdf(c_curr,c_mean',c_var); % prior probability density of current template sample
        gam_curr = interp1(tG,d_curr,t,'linear');
        SSE_curr = SSE_q(gam_curr, t, qs, q_mu_curr);
        logl_cur = f_logl(gam_curr, t, qs, q_mu_curr, sigma_cur, SSE_curr);
        
        %%%%% template function MH update %%%%%
        for ic = 1:K_mh_c
            c_prop = mvnrnd(c_curr, stepvar_chain); % Propose a candidate basis coefficient of template
            q_star_prop = c_prop*q_basis';   
            SSE_prop = SSE_q(gam_curr, t, qs, q_star_prop); % Derive the SSE using the candidate template
            logl_prop = f_logl(gam_curr, t, qs, q_star_prop, sigma_cur, SSE_prop);
            prior_c_prop = logmvnpdf(c_prop,c_mean',c_var); % prior probability density of candidate template sample
            ratio_c = logl_prop+prior_c_prop-(logl_cur+prior_c_curr);     
            if log(rand) < ratio_c
                q_mu_curr = q_star_prop;
                c_curr = c_prop;
                SSE_curr = SSE_prop;
                prior_c_curr = prior_c_prop;
                logl_cur = logl_prop;
            end         
        end
        l_pi_cur = zeros(N_func,1);
        
        %%%%% gamma block %%%%%
        for ii = 1:N_func
            kappa_ii = kappa_prop(ii);
            increments_cur = mat*d_curr(:,ii);
            l_pi_cur(ii) = sum(log(increments_cur + eps).*(kappa_prior*mat*tG));
            q_ii = qs(:,ii);

            accept_gam = 0;
            for id = 1:K_mh_d
                if((mod(id,10) == 0) && (id < 2*K_mh_d/3))
                    accept_prop = accept_gam/id;
                    if((accept_prop)>.44 ||(accept_prop)<.34)
                        kappa_ii = max(10,kappa_ii*.39/(accept_prop + 1e-10));
                        kappa_ii = min(kappa_ii, 300);
                    end
                end
            % propose candidate gamma    
            increments_star = gamrnd(mat*tG*kappa_ii,1);
            increments_star = increments_star./sum(increments_star);
            gamma_star = [0;cumsum(increments_star)];
            d_prop = interp1(tG,d_curr(:,ii),gamma_star,'linear');
            d_prop(end) = 1;
            increments_can = mat*d_prop;
            % Evaluate acceptance probability
            % evaluate priors
            l_pi_can = sum(log(increments_can + eps).*(kappa_prior*mat*tG));
            
            % evaluate likelihoods
            gam_prop = interp1(tG,d_prop,t,'linear');
            gam_prop(end) = 1;
            
            SSE_prop = SSE_q(gam_prop,t,q_ii,q_mu_curr);
            
            star_inverse = interp1(gamma_star,tG,tG,'linear');
            % Evaluate proposal densities
            star_inverse(end) = 1;
            increments_star_inverse = mat*star_inverse;
            increments_star = mat*gamma_star;
            l_g_cur_given_can = sum(log(increments_star_inverse + eps).*(kappa_ii*mat*tG));
            l_g_can_given_cur = sum(log(increments_star + eps).*(kappa_ii*mat*tG));
            % Evaluate Jacobian 
            grad_gamma_cur_inv = gradient(interp1(tG,tG,d_curr(:,ii)),tG);
            grad_gamma_can_inv = gradient(interp1(tG,tG,d_prop),tG);
            l_J = sum(log(interp1(tG,grad_gamma_can_inv,d_curr(2:(end-1),ii))))-sum(log(interp1(tG,grad_gamma_cur_inv,d_prop(2:(end-1)))));
            % Evaluate Acceptance Probability
            ratio_g = (SSE_curr(ii)-SSE_prop)/(2*sigma_cur) + (l_pi_can+l_g_cur_given_can)-(l_pi_cur(ii)+l_g_can_given_cur) + l_J;
            
            if log(rand) <= ratio_g
                d_curr(:,ii) = d_prop;
                SSE_curr(ii) = SSE_prop;
                l_pi_cur(ii) = l_pi_can;
                accept_gam = accept_gam + 1;
            end
            end
        end

        %%%% Centering Step %%%%%
        gamma_mean = SqrtMeanInverse(d_curr');
        l_pi_cur_update = zeros(N_func,1);
        for ii = 1:N_func
            d_curr(:,ii) = interp1(tG, d_curr(:,ii)', (tG(end)-tG(1)).*gamma_mean + tG(1));
            increments_cur = mat*d_curr(:,ii);
            l_pi_cur_update(ii) = sum(log(increments_cur + eps).*(kappa_prior*mat*tG));
        end
        gamma_mean = interp1(tG,gamma_mean,t);
        gam_curr = interp1(tG,d_curr,t,'linear');
        q_mu_curr = warp_q_gamma(q_mu_curr,gamma_mean,t);
        c_curr = q_mu_curr*q_basis/(q_basis'*q_basis);
        prior_q_cur_update = logmvnpdf(c_curr,c_mean',c_var);
        q_mu_curr = c_curr*q_basis';     
        SSE_curr = SSE_q(gam_curr, t, qs, q_mu_curr);

        %%%% Update weights after centering%%%%%
        w_samples(i) = w_samples(i)*exp(prior_q_cur_update-prior_c_curr+sum(l_pi_cur_update)-sum(l_pi_cur));
                
        %%%%% sigma update %%%%%
        alpha_conditional = alpha_sigma + length(t)*N_func/2;
        beta_conditional = beta_sigma + sum(SSE_curr)/2;
        invS = gamrnd(alpha_conditional,1/beta_conditional);
        sigma_cur = 1/invS;

        gam_curr = interp1(tG,d_curr,t,'linear');
        logl_cur = f_logl(gam_curr, t, qs, q_mu_curr, sigma_cur, SSE_curr);
        logpost_update(i) = prior_q_cur_update + logl_cur + sum(l_pi_cur_update) + log_invgamma_pdf(sigma_cur,alpha_sigma,beta_sigma);
        sigma_samples(i) = sigma_cur;
        c_samples(i,:) = c_curr;
        d_samples(i,:,:) = d_curr;
    end
    
    w_samples = w_samples/sum(w_samples);
    ESS(2) = 1/sum(w_samples.^2);

    w_update = w_samples;
    c_update = c_samples;
    d_update = d_samples;
    sigma_squared_update = sigma_samples;
    SMC_time = toc;
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

function [ YI ] = lsq_lut_piecewise( x, y, XI )
% LSQ_LUT_PIECEWISE Piecewise linear interpolation for 1-D interpolation (table lookup)
%   YI = lsq_lut_piecewise( x, y, XI ) obtain optimal (least-square sense)
%   vector to be used with linear interpolation routine.
%   The target is finding Y given X the minimization of function 
%           f = |y-interp1(XI,YI,x)|^2
%   
%   INPUT
%       x measured data vector
%       y measured data vector
%       XI break points of 1-D table
%
%   OUTPUT
%       YI interpolation points of 1-D table
%           y = interp1(XI,YI,x)
%

if size(x,2) ~= 1 && size(x,1) == 1
    x = x'; 
end

if size(y,2) ~= 1 && size(y,1) == 1
    y = y';
end

% matrix defined by x measurements
A = sparse([]); 

% vector for y measurements
Y = []; 

for j=2:length(XI)
    
    % get index of points in bin [XI(j-1) XI(j)]
    ix = x>=XI(j-1) & x<XI(j);

    % get x and y data subset
    x_ = x(ix);
    y_ = y(ix);
    
    % create temporary matrix to be added to A
    tmp = [(( -x_+XI(j-1) ) / ( XI(j)-XI(j-1) ) + 1) (( x_-XI(j-1) ) / ( XI(j)-XI(j-1) ))];
    
    % build matrix of measurement with constraints
    [m1,n1]=size(A);
    [m2,n2]=size(tmp);
    A = [[A zeros(m1,n2-1)];[zeros(m2,n1-1) tmp]];
    
    % concatenate y measurements of bin
    Y = [Y; y_];
end

% obtain least-squares Y estimation
YI=A\Y;
end
