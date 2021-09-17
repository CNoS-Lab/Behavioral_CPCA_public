function p = get_p_val(null_dist, val)

null_dist = null_dist(:);
for ii = 1:size(val,1)
    for jj = 1:size(val,2)
        p(ii,jj) = sum(null_dist > val(ii,jj))./length(null_dist);
    end
end

end