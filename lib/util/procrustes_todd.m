%------------------------------------------------------------------

function [B, T] = procrustes_todd(A, Target, type)
% Procrustes rotation of loadings matrix from factor analysis.

[d, m] = size(A);

if nargin < 2 || isempty(Target)
    error('A target matrix must be specified for procrustes rotation.');
elseif any(size(Target) ~= [d m])
    error('Incorrect size for procrustes rotation target matrix.');
end
if nargin < 3 || isempty(type)
    type = 'oblique';
else
    typeNames = {'oblique','orthogonal'};
    i = strmatch(lower(type), typeNames);
    if length(i) > 1
        error(sprintf('Ambiguous ''TypeProcr'' parameter value:  %s.', type));
    elseif isempty(i)
        error(sprintf('Unknown ''TypeProcr'' parameter value:  %s.', type));
    end
    type = typeNames{i};
end

% Orthogonal rotation to target
switch type
case 'orthogonal'
    [L, D, M] = svd(Target' * A);
    T = M * L';
    
% Oblique rotation to target
case 'oblique'
    % LS, then normalize
    T = A \ Target;
    T = T * diag(sqrt(diag((T'*T)\eye(m)))); % normalize inv(T)
end
B = A * T;
