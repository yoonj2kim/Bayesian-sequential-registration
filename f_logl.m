function out = f_logl(gamma, t, qs, q_mu, sigma_squared, SSEg)
    % sum of log observation likelihood
    if (SSEg == 0)
        SSEg = SSE_q(gamma, t, qs, q_mu);
    end
    n = length(t);
    out = n * log(1/sqrt(2*pi)) - n * log(sqrt(sigma_squared)) - (SSEg ./ (2 * sigma_squared));
    out = sum(out);
end
