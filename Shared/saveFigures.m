function [] = saveFigures(figures, savePath)
    
    cellfun(@(fig) saveas(fig.handle, [savePath, filesep, 'plot_', fig.name, '.fig']), figures);
    cellfun(@(fig) saveas(fig.handle, [savePath, filesep, 'plot_', fig.name, '.svg']), figures);
    
end