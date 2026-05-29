function plotScalingFactors(modParam_scaling, modParam_names, nModParams, savePath)    
%  plot scaling factor distributions

    fig = {figure('Name','scalingFactors')}; set(gcf, 'Position', [100, 100, 1200, 1000]);
    for iModParam = 1:nModParams
        subplot(4, ceil(nModParams/4), iModParam);
        iScalingDistribution = modParam_scaling(:,iModParam);
        plotDist = iScalingDistribution>0;
        if isempty(find(iScalingDistribution(plotDist)>1))
            % Define bin edges to ensure a single vertical line at 1
            binEdges = [0.99, 1.01]; % Bins that capture the value 1
            histogram(iScalingDistribution(plotDist), 'BinEdges', binEdges,'FaceColor', 'blue');
            xlim([0.5,1.5])
        else
            histogram(iScalingDistribution(plotDist), 'FaceColor', 'blue');
        end
        title(modParam_names{iModParam});
        grid on;
        if iModParam > nModParams-4
            xlabel('scaling factor value')
        end
        if iModParam-1 == 0 || mod(iModParam-1,4) == 0
            ylabel('number of cells')
        end
    end
    sgtitle('Scaling Factor Distributions', 'Interpreter', 'none');
    figname = strcat(savePath,filesep,'Distributions');
    % savefig(figname);
end