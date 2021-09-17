function out_thresh = est_corr_thresholds(n_dim,n_bootstrap,p_val)

r_val = NaN(n_bootstrap,1);

for ii = 1:n_bootstrap
    r_val(ii) = corr(randn(n_dim,1),randn(n_dim,1));
end

out_thresh = quantile(abs(r_val), 1-p_val);

end