function log_pdf = log_invgamma_pdf(x, alpha, beta)
    log_pdf = alpha * log(beta) - gammaln ( alpha ) - (alpha + 1) .* log (x) - beta ./ x;
end
