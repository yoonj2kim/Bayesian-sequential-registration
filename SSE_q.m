function out = SSE_q(gamma,t,qs,q_mu)
    % sum of squared difference between q_i and (q_\mu,\gamma_i^{-1})
    % across t
    out = zeros(1,size(gamma,2));
    for ii = 1:size(qs,2)
        warped_q_mean = warp_q_gamma(q_mu,invertGamma(gamma(:,ii)'),t);
        vec = (qs(:,ii)-warped_q_mean).^2;
        out(ii) = sum(vec);
    end
end