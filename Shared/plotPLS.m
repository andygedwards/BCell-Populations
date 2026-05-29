function R2_output = plotPLS(PLS_output, X_vals_raw, Y_vals_raw, Y_vals_log, X_names, Y_names, X_model, Y_model, analysis_mode, analysis_direction, analysis_processing, savePath)

    %% Function Input Validation
    if ~any(strcmpi(analysis_mode, {'sensitivity','translation'}))
        error('Analysis mode must be "sensitivity" or "translation".');
    end
    if ~any(strcmpi(analysis_direction, {'forward','reverse'}))
        error('Analysis direction must be "forward" or "reverse".');
    end
    if ~any(strcmpi(analysis_processing, {'raw','normalized'}))
        error('Analysis processing must be "raw" or "normalized".');
    end
    
    %% Print Program Status to Command Window
    if strcmpi(analysis_mode, 'sensitivity')
        fprintf("Plotting %s %s %s analysis results... ", X_model, lower(analysis_direction), lower(analysis_mode));
    elseif strcmpi(analysis_mode, 'translation')
        fprintf("Plotting %s-%s %s analysis results... ", X_model, Y_model, lower(analysis_mode));
    end
    runTime = tic;
    
    %% Parent Folder
    
    if strcmpi(analysis_mode, 'sensitivity')
        figureFolder_parent = [savePath, filesep, 'Sensitivity Analyses'];
    elseif strcmpi(analysis_mode, 'translation')
        figureFolder_parent = [savePath, filesep, 'Translation Analyses'];
    end
    if exist(figureFolder_parent, 'dir') ~= 7
        mkdir(figureFolder_parent);
    end

    %% Parse NIPALS Output

    % T = PLS_output{1};
    % P = PLS_output{2};
    % W = PLS_output{3};
    % W_star = PLS_output{4};
    % U = PLS_output{5};
    % B = PLS_output{6};
    % C = PLS_output{7};
    bPLS = real(PLS_output{8});
    % bPLS_star = PLS_output{9};
    % X_hat = PLS_output{10};
    Y_hat = PLS_output{11};
    % X_R2 = PLS_output{12};
    % Y_R2 = PLS_output{13}; 

    %% Define Plotting Settings

    pointSize = 7;
    X_nVars = size(X_vals_raw, 2); 
    Y_nVars = size(Y_vals_raw, 2); 
    nSubplots = 6; % number of subplots per figure
    Y_nFigures = ceil(Y_nVars / nSubplots); 
    X_nFigures = ceil(X_nVars / nSubplots);

    %% Calculate Statistics

    SSYT = real(sum((Y_vals_log - mean(Y_vals_log)) .^ 2)); % actual measured outputs, nonregressed (predicted/true)
    SSYR = real(sum((Y_hat - mean(Y_vals_log)) .^ 2)); % Predicted outcome for each feature of interest (residual)
    R2_each = SSYR ./ SSYT; 
    R2_avg = mean(R2_each);
    R2_std = std(R2_each);
    R2_output = {R2_avg, R2_std};
    Y_hat_real = real(exp(Y_hat));

    X_vals_avg = mean(X_vals_raw);
    X_vals_dev = std(X_vals_raw);
    Y_vals_avg = mean(Y_vals_raw);
    Y_vals_dev = std(Y_vals_raw);

    if strcmpi(analysis_mode, 'sensitivity') && strcmpi(analysis_direction, 'forward')
        %% Create Output Folder for Figures
        figureFolder_stats = [savePath, filesep, 'Population Stats', filesep, X_model, '_Population Stats'];
        if exist(figureFolder_stats, 'dir') ~= 7
            mkdir(figureFolder_stats);
        end
        figures = {};
        
        %% Plot Input Parameter (X) and Output Feature (Y) Distributions
        
        % Parameter Distribution
        iVar = 1;
        for iFigure = 1:X_nFigures
%             figureFile = [figureFolder, 'plot_', ];
%             if exist(figureFolder_general, 'dir') ~= 7
%                 mkdir(figureFolder_translation);
%             end
            figures = [figures, {newFigure([X_model, '_Parameter Distribution-', num2str(iFigure)])}];
            set(gcf, 'Position', [50,100,1500,750]);
            for iSubplot = 1:nSubplots
                if iVar <= X_nVars
                    subplot(2,3,iSubplot);
                    rows = find(X_vals_raw(:,iVar)~=0);
                    X_vals_avg(iVar) = mean(X_vals_raw(rows,iVar));
                    X_vals_dev(iVar) = std(X_vals_raw(rows,iVar));
                    histogram(X_vals_raw(rows,iVar), 'FaceColor', 'blue') ;
                    set(gca, 'box', 'off', 'tickdir',' out', 'fontsize', 10);
                    xlabel([X_names{iVar}], 'Interpreter', 'none');
                    title(['Mean = ', num2str(X_vals_avg(iVar),3), '; Std = ', num2str(X_vals_dev(iVar),3)], 'Interpreter', 'none')
                    grid on;
                    iVar = iVar + 1;
                end
            end
            sgtitle([X_model, ': Parameter Distributions']);
        end

        % Feature Distribution
        iVar = 1;
        for iFigure = 1:Y_nFigures
            figures = [figures, {newFigure([X_model, '_Feature Distribution-', num2str(iFigure)])}];
            set(gcf, 'Position', [50,100,1500,750]);
            for iSubplot = 1:nSubplots
                if iVar <= Y_nVars
                    subplot(2,3,iSubplot);
                    histogram(Y_vals_raw(:,iVar), 'FaceColor', 'blue') ;
                    set(gca, 'box', 'off', 'tickdir',' out', 'fontsize', 10);
                    xlabel([Y_names{iVar}], 'Interpreter', 'none');
                    title(['Mean = ', num2str(Y_vals_avg(iVar),3), '; Std = ', num2str(Y_vals_dev(iVar),3)], 'Interpreter', 'none');
                    grid on;
                    iVar = iVar + 1;
                end
            end
            sgtitle([Y_model, ': Feature Distributions']);
        end
        
        %% Feature-Feature (Y-Y) Correlations 
    
        iSubplot = 1;
        figures = [figures, {newFigure([X_model, '_Feature Correlations'])}];
        for iVar = 1:Y_nVars
            for jVar = 1:Y_nVars
                subplot(Y_nVars, Y_nVars, iSubplot);
                set(gca,'box', 'off', 'tickdir', 'out', 'fontsize', 10);
                iSubplot = iSubplot + 1;
                scatter(Y_vals_raw(:,jVar), Y_vals_raw(:,iVar), pointSize/2, 'filled', 'blue');
                grid on;
                if iVar == Y_nVars
                    xlabel(Y_names{jVar}, 'Interpreter', 'none');
                end
                if jVar == 1
                    ylabel(Y_names{iVar}, 'Interpreter', 'none');
                end
                R = corrcoef(Y_vals_raw(:,jVar), Y_vals_raw(:,iVar));
                title(['R^2 = ', num2str(R(1,2)^2,3)]);
            end
        end
        sgtitle([Y_model, ': Feature Correlations']);
        
        saveFigures(figures, figureFolder_stats);
    end
    
    %% Create Output Folder for Figures
    
    if strcmpi(analysis_mode, 'sensitivity')
        figureFolder_child = [figureFolder_parent, filesep, X_model, '_Sensitivity_', analysis_direction];
    elseif strcmpi(analysis_mode, 'translation')
        figureFolder_child = [figureFolder_parent, filesep, X_model, '-', Y_model, '_Translation'];
    end
    if exist(figureFolder_child, 'dir') ~= 7
        mkdir(figureFolder_child);
    end
    figures = {};

    %% Plot PLS Predictions

    iVar = 1;
    for iFigure = 1:Y_nFigures
        if strcmpi(analysis_mode, 'sensitivity')
            figures = [figures, {newFigure([X_model, '_', analysis_direction, '_Sensitivity Prediction-', num2str(iFigure)])}];
        elseif strcmpi(analysis_mode, 'translation')
            figures = [figures, {newFigure([X_model, '-', Y_model, '_Translation Prediction-', num2str(iFigure)])}];
        end
        set(gcf, 'Position', [50,100,1500,750]);
        for iSubplot = 1:nSubplots
            if iVar <= Y_nVars 
                subplot(2,3,iSubplot);
                hold on;
                
                % Plot data points
                scatter(Y_vals_raw(:,iVar), Y_hat_real(:,iVar), pointSize, 'filled', 'blue');
                if strcmpi(analysis_mode, 'sensitivity')
                    xlabel('Observed');
                    ylabel('Predicted');
                elseif strcmpi(analysis_mode, 'translation')
                    xlabel([Y_model, ' - Observed']);
                    ylabel([Y_model, ' - Predicted']); 
                end
                
                title([Y_names{iVar}, ' (R2 = ', num2str(R2_each(iVar), 4), ')'], 'Interpreter', 'none');
                set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 10);
                grid on;
                
                % Plot identity line
                ylim_ind = get(gca, 'ylim');
                xlim_ind = get(gca, 'xlim');
                minpoint = min([ylim_ind(1),xlim_ind(1)]);
                maxpoint = max([ylim_ind(2),xlim_ind(2)]);
                plot([minpoint, maxpoint], [minpoint,maxpoint], '--k', 'LineWidth', 2);
                xlim([minpoint,maxpoint]);
                ylim([minpoint,maxpoint]);
                try
                    axis equal;
                catch
                    warning("Unable to create equal plot axis scaling.");
                end
                iVar = iVar + 1;
                hold off;  
            end
        end
        if strcmpi(analysis_mode, 'sensitivity')
            sgtitle([X_model, ': ', analysis_direction, ' Sensitivity Prediction']);
        elseif strcmpi(analysis_mode, 'translation')
            sgtitle([X_model, '\rightarrow', Y_model, ': Translation Prediction']);
        end
    end

    %% Plot PLS Error

    if strcmpi(analysis_mode, 'translation')
        err = (Y_hat_real - Y_vals_raw) ./ Y_vals_raw;
        errAvg = mean(abs(err));
        errStd = std(abs(err));
        iVar = 1;
        for iFigure = 1:Y_nFigures
            figures = [figures, {newFigure([X_model, '-', Y_model,  '_', '_Translation Error (scat)-', num2str(iFigure)])}];
            for iSubplot = 1:nSubplots
                if iVar <= Y_nVars 
                    subplot(2,3,iSubplot);
                    scatter(X_vals_raw(:,iVar), err(:,iVar), pointSize, 'filled', 'blue'); 
                    xlabel([X_model, ' - Observed']);
                    ylabel([Y_model ,' - Prediction Error (%)']);
                    title([Y_names{iVar}, ' (avg = ', num2str(errAvg(iVar),2), '; std = ', num2str(errStd(iVar),2), ')'], 'Interpreter', 'none');
                    set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 10);
                    grid on;
                    axis normal;
                    iVar = iVar + 1;
                end
            end
            sgtitle([X_model, '\rightarrow', Y_model, ': Translation Error']);
        end

        figures = [figures, {newFigure([X_model, '-', Y_model,  '_', '_PLS Error (bar)'])}];
        hold on;
        bar(errAvg, 'FaceColor', 'blue');
        errorbar(errAvg,errStd,'k','linestyle','none','LineWidth',1);
        set(gca, 'XTick', 1:Y_nVars)
        set(gca, 'XTickLabel', Y_names, 'TickLabelInterpreter', 'none')
        ylabel("Mean Absolute Prediction Error (%)");
        ylim([0,inf]);
        title([X_model, '\rightarrow', Y_model,': Translation Error']);
        grid on;
        hold off;
    end

    %% Plot PLS R-Squared Values 
    if strcmpi(analysis_mode, 'sensitivity')
        figures = [figures, {newFigure([X_model, '_', analysis_direction, '_Sensitivity R2'])}];
    elseif strcmpi(analysis_mode, 'translation')
        figures = [figures, {newFigure([X_model, '-', Y_model, '_Translation R2'])}];
    end
    bar(R2_each, 'FaceColor', 'blue')
    set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 12)
    set(gca, 'XTick', 1:Y_nVars)
    set(gca, 'XTickLabel', Y_names, 'TickLabelInterpreter', 'none')
    set(gca, 'XLim', [0,Y_nVars+1])
    ylim([0, 1])
    grid on;
    if strcmpi(analysis_mode, 'sensitivity')
        title([X_model, ': ', analysis_direction, ' PLS R^2 values (mean = ', num2str(R2_avg,3), ')']);
    elseif strcmpi(analysis_mode, 'translation')
        title([X_model, '\rightarrow', Y_model, ': PLS R^2 values (mean = ', num2str(R2_avg,3), ')']);
    end
    
    %% Plot bPLS Matrix Coefficients 
    
    % Grouped By Y
    iVar = 1;
    for iFigure = 1:Y_nFigures
        if strcmpi(analysis_mode, 'sensitivity')
            figures = [figures, {newFigure([X_model, '_', analysis_direction, '_Sensitivity Coefficients -Y', num2str(iFigure)])}];
        elseif strcmpi(analysis_mode, 'translation')
            figures = [figures, {newFigure([X_model, '-', Y_model, '_Translation Coefficients -Y', num2str(iFigure)])}];
        end
        set(gcf, 'Position', [50,100,1500,750]);
        for iSubplot = 1:nSubplots
            if iVar <= Y_nVars
                subplot(2, 3, iSubplot);
                bar(bPLS(:,iVar), 'FaceColor', 'blue');
                title(['(Y) ', Y_names{iVar}], 'Interpreter', 'none');
                set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 10);
                set(gca, 'XTick', 1:X_nVars);
                set(gca, 'XLim', [0 X_nVars+1]);
                set(gca, 'YLim', [-1,1]);
                set(gca, 'XTickLabel', X_names, 'TickLabelInterpreter', 'none');
                grid on;
                iVar = iVar + 1;
            end
        end
        if strcmpi(analysis_mode, 'sensitivity')
            sgtitle([X_model, ': ', analysis_direction, ' PLS Coefficients']);
        elseif strcmpi(analysis_mode, 'translation')
            sgtitle([X_model, '\rightarrow', Y_model, ': PLS Coefficients']);
        end
    end
    
    % Grouped By X
    iVar = 1;
    for iFigure = 1:X_nFigures
        if strcmpi(analysis_mode, 'sensitivity')
            figures = [figures, {newFigure([X_model, '_', analysis_direction, '_Sensitivity Coefficients -X', num2str(iFigure)])}];
        elseif strcmpi(analysis_mode, 'translation')
            figures = [figures, {newFigure([X_model, '-', Y_model, '_Translation Coefficients -X', num2str(iFigure)])}];
        end
        set(gcf, 'Position', [50,100,1500,750]);
        for iSubplot = 1:nSubplots
            if iVar <= X_nVars
                subplot(2, 3, iSubplot);
                bar(bPLS(iVar,:), 'FaceColor', 'blue'); 
                title(['(X) ', X_names{iVar}], 'Interpreter', 'none');
                set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 10);
                set(gca, 'XTick', 1:Y_nVars);
                set(gca, 'XLim', [0,Y_nVars+1]);
                set(gca, 'YLim', [-1,1]);
                set(gca, 'XTickLabel', Y_names, 'TickLabelInterpreter', 'none');
                grid on;
                iVar = iVar + 1;
            end
        end
        if strcmpi(analysis_mode, 'sensitivity')
            sgtitle([X_model, ': ', analysis_direction, ' PLS Coefficients']);
        elseif strcmpi(analysis_mode, 'translation')
            sgtitle([X_model, '\rightarrow', Y_model, ': PLS Coefficients']);
        end
    end
    
    saveFigures(figures, figureFolder_child);
    fprintf("%.3fs\n", toc(runTime));

end