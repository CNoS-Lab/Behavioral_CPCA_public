function [Z, G] = removeNaN(Z, G)

% This function removes rows from Z and G that contain NaNs

temp = [Z, G];
idx = isnan(mean(temp,2));
Z(idx,:) = [];
G(idx,:) = [];

end