function [h] = plotLoadings(X,xTicks,yTicks)

h = heatmap(X);
h.XDisplayLabels = xTicks;
h.YDisplayLabels = yTicks;

end