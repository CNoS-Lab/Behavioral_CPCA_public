function [ax] = plotLoadings2(X,xTicks,yTicks,thresh)

h = imagesc(X);
xticklabels(xTicks);
yticklabels(yTicks);

ax = gca;
ax.XTick = 1:size(X,2);
ax.XTickLabel = xTicks;
ax.YTick = 1:size(X,1);
ax.YTickLabel = yTicks;

[rows,cols] = size(X);

for i = 1:rows
    for j = 1:cols
        if X(i,j) >= thresh
            textHandles(j,i) = text(j,i,[num2str(round(X(i,j),2)),'*'],...
                'horizontalAlignment','center');
        else
            textHandles(j,i) = text(j,i,num2str(round(X(i,j),2)),...
                'horizontalAlignment','center');
        end
    end
end

ax = gca;

end