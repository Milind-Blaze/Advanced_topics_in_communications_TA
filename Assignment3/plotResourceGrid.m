%{ 
Function to plot the generated resource grid

plotResourceGrid(data, titleValue, xlabelValue, ylabelValue)

Inputs:
    data (2D array): data/RG for which heatmap needs to be generated
    titleValue (string): title of the plot
    xlabelValue (string): X axis label of the figure
    ylabelValue (string): Y axis label of the figure
%}

function plotResourceGrid(data, titleValue, xlabelValue, ylabelValue)
    figure
    imagesc(data);
    colormap('jet');
    title(titleValue)
    colorbar;
    set(gca,'XTick',[1:14])
    set(gca,'YDir','normal')
    xlabel(xlabelValue);
    ylabel(ylabelValue);
    
end