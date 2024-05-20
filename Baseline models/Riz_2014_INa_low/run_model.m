
% Set up model parameters
[param, param_names] = Riz2014_init_parameters_INa_low();
% Set up initial conditions
[states, state_names] = Riz2014_init_states_INa_low();

param(find(strcmp(param_names, 'glycolysis'))) = 1;

% Low glucose run
param(find(strcmp(param_names, 'G'))) = 2;

Tstop = 6000000;  % simulation time
options = [];
[T_init, Y_init] = ode15s(@Riz2014_rhs_INa_low, [0, Tstop], states, options, param);

Tstop = 10000000;  % simulation time
[T_low, Y_low] = ode15s(@Riz2014_rhs_INa_low, [0, Tstop], Y_init(end,:), options, param);

% High glucose run
param(find(strcmp(param_names, 'G'))) = 10;

[T_high, Y_high] = ode15s(@Riz2014_rhs_INa_low, [0, Tstop], states, options, param);

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
title('2 mM Glucose - glycolysis decoupled')

subplot(2,1,2), plot(T_high/1000, Y_high(:,find(strcmp(state_names, 'V'))), 'linewidth', 2)
set(gca, 'fontsize', 14)
xlabel('t (sec)')
ylabel('mV')
ylim([-80,0])
title('10 mM Glucose - glycolysis decoupled')

