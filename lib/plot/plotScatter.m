function [] = plotScatter(x, y, Str, labelsStr, titleStr)

scatter(x,y,80,Str,'filled');
dx = 0.02; dy = 0.02;

if nargin > 3
    suptitle(titleStr);
    text(x+dx, y+dy, labelsStr);
elseif nargin > 2
    text(x+dx, y+dy, labelsStr);
end

xlabel('Liberal SgD/ Perceptual distortions');
ylabel('Distraction on DL/ Sensory overload');
axis square; 

end