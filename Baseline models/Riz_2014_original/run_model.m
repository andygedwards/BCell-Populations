% test
% Set up model parameters
[param, param_names] = Riz2014_init_parameters();
% Set up initial conditions
[states, state_names] = Riz2014_init_states();

param(find(strcmp(param_names, 'glycolysis'))) = 1; %Param 41 set to 1

% Low glucose run
param(find(strcmp(param_names, 'G'))) = 2; %Param 51 set to 2

Tstop = 6000000;  % simulation time %ms
options = [];
[T_init, Y_init] = ode15s(@Riz2014_rhs, [0, Tstop], states, options, param); %Run low glu to SS

Tstop = 10000000;  % simulation time
[T_low, Y_low] = ode15s(@Riz2014_rhs, [0, Tstop], Y_init(end,:), options, param); %Use value at t126 

% High glucose run
param(find(strcmp(param_names, 'G'))) = 10;

%Run for high glucose
[T_high, Y_high] = ode15s(@Riz2014_rhs, [0, Tstop], states, options, param); %Use SS value for high glucose?

% Set up figure
 figure('Units','centimeters', 'Position', [2 10 15 12], ...
    'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'PaperSize', [15, 12]);

% Plot the action potential
subplot(2,1,1), plot(T_low/1000, Y_low(:,find(strcmp(state_names, 'V'))), 'linewidth', 2)
set(gca, 'fontsize', 14)
xlabel('t (sec)')
ylabel('mV')
ylim([-80,0])
title('2 mM Glucose')

subplot(2,1,2), plot(T_high/1000, Y_high(:,find(strcmp(state_names, 'V'))), 'linewidth', 2)
set(gca, 'fontsize', 14)
xlabel('t (sec)')
ylabel('mV')
ylim([-80,0])
title('10 mM Glucose')

%% Plot currents vs voltage
current_matrix = zeros(length(T_high), 13);
ylabels_sub = {'mV', 'INa', 'IBK', 'ICaL', 'ICaPQ', 'ICaT', 'IGABAR', 'IHERG', 'IKATP', 'IKv', 'ISK', 'Ileak', 'Istim', 'Ivclamp'};

for i = 1: length(T_high)
    [~, currents2] = Riz2014_rhs(T_high(i), Y_high(i,:), param);
    current_matrix(i, :) = currents2;

end

% Why HERG channel conductance is set to 0 ?

figure(3); hold on
index_plot = find(T_high./1000 > 9860 & T_high./1000 < 9900 );
subplot(7,2,1), 
plot(T_high(index_plot)/1000, Y_high(index_plot,14), 'linewidth', 1)
set( ...
    gca, 'fontsize', 14)
xlabel('t (sec)'); ylabel(ylabels_sub{1});

for j = 2:14
    subplot(7,2,j), 
    plot(T_high(index_plot)/1000, current_matrix((index_plot), j-1), 'linewidth', 1)
    set(gca, 'fontsize', 14)
    xlabel('t (sec)'); ylabel(ylabels_sub{j}); 
end

figure(4); hold on
for j = 1:14
    subplot(7,2,j), 
    plot(T_high(index_plot)/1000, Y_high(index_plot,j), 'linewidth', 1)
    set(gca, 'fontsize', 14)
    xlabel('t (sec)'); ylabel(state_names{j}); 
end




