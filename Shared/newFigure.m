function [figures] = newFigure(saveName)

    figures.handle = figure('visible','off','NumberTitle','off');
    set(gcf,'Visible','off','CreateFcn','set(gcf,''Visible'',''on'')') %save, don't show.
    figures.name = saveName;

end