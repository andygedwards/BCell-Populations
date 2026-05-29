function [] = plot_refined(Na,Combined,Na_cost,Combined_cost,variable)
    % Create figure with adjusted width
    figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.3]);
    
    % First subplot: Bar plot with error bars and group names
    subplot(1,3,1);
    data = [Na; Combined]'; % Combine data into a matrix
    means = mean(data);              % Compute means
    std_devs = std(data);            % Compute standard deviations
    bar(means);                      % Bar plot
    hold on;
    errorbar(1:length(means), means, std_devs, 'k', 'linestyle', 'none'); % Add error bars
    hold off;
    xticks(1:length(means));          % Set x-tick positions
    xticklabels({'Na', 'Combined'}); % Set group names
    ylabel(['Mean ',variable]);             % Label for y-axis
    
    % Second subplot: Scatter plot for Combined vs. Combined_cost
    subplot(1,3,2);
    scatter(Combined, Combined_cost);
    xlabel(['Combined ',variable]);         % Label for x-axis
    ylabel('Combined Cost');         % Label for y-axis
    
    % Third subplot: Scatter plot for Na vs. Na_cost
    subplot(1,3,3);
    scatter(Na, Na_cost);
    xlabel(['Na ', variable]);               % Label for x-axis
    ylabel('Na Cost');               % Label for y-axis
end