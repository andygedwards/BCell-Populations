function [early_exocytosis, late_exocytosis, total_exocytosis] = Exocytosis(parameters, parameter_names, states, state_names, rhs_model)
    % Camunas-Soler protocols
    parameters(find(strcmp(parameter_names, 'g_KATP0'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_KATP_hat'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_SK'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_BK'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_Kv'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_HERG'))) = 0;
    parameters(find(strcmp(parameter_names, 'V_K'))) = -92;
    
    % Multi-pulse current/capacitance
    v_rest = -70;
    v_pulse = 0;
    prepulse_length = 1000;
    end_length = 1000;
    step_length = 500;
    cycle_length = 1000;
    num_pulses = 10;
    v_clamp_amps(1:2:2*num_pulses+1) = v_rest; v_clamp_amps(2:2:num_pulses*2+1) = v_pulse;
    v_clamp_times(1) = prepulse_length; v_clamp_times(2:num_pulses*2+1) = (prepulse_length+step_length:step_length:prepulse_length+num_pulses*2*step_length); v_clamp_times(end) = v_clamp_times(num_pulses*2+1)+end_length; 
    v_clamp_times = [5000 v_clamp_times(2:end) + 4000];
    R_clamp = 0.001; % (GOhm)
    
    
    current = 'exocytosis';
    plotting = 'true';
    options = odeset('MaxStep',1);
    
    init_tspan = [0,10000];
    [T_init, Y_init] = ode15s(rhs_model, init_tspan, states, options, parameters, v_rest, init_tspan(end), R_clamp);
    
    tspan = [0,v_clamp_times(end)];
    [T, Y] = ode15s(rhs_model, tspan, Y_init(end,:), options, parameters, v_clamp_amps, v_clamp_times, R_clamp);
    
    try
        outs = analyze_v_clamp(Y,T,parameters,state_names, cycle_length, prepulse_length, end_length, rhs_model, v_clamp_amps, v_clamp_times, R_clamp, current, plotting);
        early_exocytosis = outs{1};
        late_exocytosis = outs{2};
        total_exocytosis = outs{3};
    catch ME
        reason = {ME.identifier};
        early_exocytosis = NaN;
        late_exocytosis = NaN;
        total_exocytosis = NaN;
    end
end