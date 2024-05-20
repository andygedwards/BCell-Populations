function plotScalingFactors(modParam_scaling, modParam_names, nModParams, savePath)    
%  plot scaling factor distributions

    fig = {figure('Name','scalingFactors')};
    for iModParam = 1:nModParams
        subplot(4, ceil(nModParams/4), iModParam);
        iScalingDistribution = modParam_scaling(:,iModParam);
        plotDist = iScalingDistribution>0;
        histogram(iScalingDistribution(plotDist), 'FaceColor', 'blue');
        title(sprintf("%s: Mean = %.3f, Stdev = %.3f", modParam_names{iModParam}, mean(iScalingDistribution(plotDist)), std(iScalingDistribution(plotDist))), 'Interpreter', 'none');
        grid on;
    end
    sgtitle('Scaling Factor Distributions', 'Interpreter', 'none');
    figname = strcat(savePath,'/','Distributions');
    savefig(figname);
end