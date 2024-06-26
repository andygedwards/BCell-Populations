function [v_half, f] = INa_inact(parameters, parameter_names, states, state_names, rhs_model)
    % Camunas-Soler protocols
    parameters(find(strcmp(parameter_names, 'g_KATP0'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_KATP_hat'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_SK'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_BK'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_Kv'))) = 0;
    parameters(find(strcmp(parameter_names, 'g_HERG'))) = 0;
    parameters(find(strcmp(parameter_names, 'V_K'))) = -92;
    
    % Na+ availability
    v_rest = (-120:10:-20);
    v_int = -80;
    v_pulse = -10;
    rest_length = 500;
    int_length = 2;
    pulse_length = 50;
    cycle_length = rest_length+int_length+pulse_length;
    prepulse_length = 0;
    end_length = 0;
    num_pulses = length(v_rest);
    current = 'I_Na';
    plotting = 'false';
    options = odeset('MaxStep',1);
    
    v_clamp_amps(1:3:3*num_pulses) = v_rest; v_clamp_amps(2:3:num_pulses*3) = v_int; v_clamp_amps(3:3:num_pulses*3) = v_pulse;
    v_clamp_times(1:3:3*num_pulses) = rest_length:cycle_length:num_pulses*cycle_length; v_clamp_times(2:3:num_pulses*3) = int_length+rest_length:cycle_length:num_pulses*cycle_length; v_clamp_times(3:3:num_pulses*3) = cycle_length:cycle_length:num_pulses*cycle_length; 
    R_clamp = 0.001; % (GOhm)
    
    init_tspan = [0,10000];
    [T_init, Y_init] = ode15s(rhs_model,init_tspan, states, options, parameters, v_int, init_tspan(end), R_clamp);
    
    tspan = [0,v_clamp_times(end)];
    [T, Y] = ode15s(rhs_model,tspan, Y_init(end,:), options, parameters, v_clamp_amps, v_clamp_times, R_clamp);
    
    try
        outs = analyze_v_clamp(Y,T,parameters,state_names, cycle_length, prepulse_length, end_length, rhs_model, v_clamp_amps, v_clamp_times, R_clamp, current, plotting);
        v_half = outs{1};
        f = outs{2};
    catch ME
        reason = {ME.identifier};
        disp(reason)
        v_half = NaN;
        f = NaN;
    end
end