function plotdist(data, bins, name, series, frac)

    if strcmp(name, 'High and Low {\it{I_{Na}}} peak')
        for i = 1: length(data)
            if i > 1
                plot_data{i} = data{i} ./ mean(vertcat(data{2:end}));
            else
                plot_data{i} = data{i} ./ mean(data{i});
            end
        end
        range_data = abs(min(plot_data{1}) - max(plot_data{1}));
    elseif strcmp(name, 'High and Low {\it{I_{Na}}} half inactivation')
        for i = 1: length(data)
            plot_data{i} = data{i};
        end 
        range_data = abs(min(plot_data{1}) - max(plot_data{1}));
    else
        for i = 1: length(data)
            plot_data{i} = data{i} ./ mean(data{i});
            range_data(i) = abs(min(plot_data{i}) - max(plot_data{i}));
        end
    end

    % Calculate the bin width for the histogram
    range = min(range_data);
    bin_width = range / bins;  % Define a reasonable bin width
    % Colormap
    colors = {[125, 62, 119] / 255 ...  % violet
            [85, 138, 163] / 255 ...    % blue
            [246, 88, 95] / 255 ...  % red
            }; 

    % Plot the series
    for i = 1:length(series)
        % Calculate histogram counts and bin edges manually
        [counts, edges] = histcounts(plot_data{i}, 'BinWidth', bin_width, 'Normalization', 'probability');

        % Compute bin centers from the edges
        bin_centers = edges(1:end-1) + diff(edges) / 2;

        % If it's the 'High I_{Na} Model' series, scale the counts
        if strcmp(series{i}, 'High $I_{Na}$ Model')
            scaling_factor = frac;  % Adjust the scaling factor as needed
            counts = counts * scaling_factor;
        elseif strcmp(series{i}, 'Low $I_{Na}$ Model')
            scaling_factor = 1-frac;  % Adjust the scaling factor as needed
            counts = counts * scaling_factor;
        end

        % Plot the scaled histogram manually using the bar function
        bar(bin_centers, counts, 'FaceColor', colors{i}, 'FaceAlpha', 0.6, 'BarWidth', 1);
        hold on;
    end
    title(name); 
    xlabel(['Normalized ', name]); ylabel('Frequency')
end