function [] = plotScree(D,titleStr)
var = D;
k = length(D);
plot(1:k,var(1:k),'-O');
xlabel('No. of Components');
ylabel('Eig Values of nth component');
if nargin<2
    title('Scree Plot');
else
    title(titleStr);
end
end