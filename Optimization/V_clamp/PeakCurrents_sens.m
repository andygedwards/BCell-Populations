function [peak_INa, v_half_act, vha_r2, peak_ICa, late_ICa] = PeakCurrents(parameters, parameter_names, states, state_names, rhs_model)
    % Camunas-Soler protocols
    parameters(find(strcmp(parameter_names, 'g_KATP0'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_KATP_hat'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_SK'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_BK'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_Kv'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_HERG'))) = 0;
    parameters(find(strcmp(parameter_names, 'V_K'))) = -92;
    
    % Peak currents
    v_rest = -70;
    v_pulse = (-60:10:30);
    rest_length = 1000;
    pulse_length = 500;
    cycle_length = rest_length+pulse_length;
    prepulse_length = 0;
    end_length = 0;
    num_pulses = 10;
    current = 'peak currents';
    plotting = 'true';
    options = odeset('MaxStep',1);

    v_clamp_amps(1:2:2*num_pulses) = v_rest; v_clamp_amps(2:2:num_pulses*2) = v_pulse;
    v_clamp_times(1:2:2*num_pulses) = rest_length:cycle_length:num_pulses*cycle_length; v_clamp_times(2:2:num_pulses*2) = cycle_length:cycle_length:num_pulses*cycle_length; 
    v_clamp_times = [5000 v_clamp_times(2:end) + 4000];
    R_clamp = 0.001; % (GOhm)
    
    init_tspan = [0,10000];
    [T_init, Y_init] = ode15s(rhs_model,init_tspan, states, options, parameters, v_rest, init_tspan(end), R_clamp);
    
    tspan = [0,v_clamp_times(end)];
    [T, Y] = ode15s(rhs_model,tspan, Y_init(end,:), options, parameters, v_clamp_amps, v_clamp_times, R_clamp);
    
    try
        outs = analyze_v_clamp(Y,T,parameters,state_names, cycle_length, prepulse_length, end_length, rhs_model, v_clamp_amps, v_clamp_times, R_clamp, current, plotting);
        peak_INa = outs{1};
        v_half_act = outs{2};
        vha_r2 = outs{3};
        peak_ICa = outs{4};
        late_ICa = outs{5};
    catch ME
        reason = {ME.identifier};
        disp(reason)
        peak_INa = NaN;
        peak_ICa = NaN;
        late_ICa = NaN;
    end
end