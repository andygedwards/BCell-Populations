function [ramp_IK] = HerringtonRamp(parameters, parameter_names, states, state_names, rhs_model)
    % Herrington protocols
%     parameters(find(strcmp(parameter_names, 'g_KATP0'))) = 0.01;
%     parameters(find(strcmp(parameter_names, 'g_KATP_hat'))) = 0.05;
%     parameters(find(strcmp(parameter_names, 'g_SK'))) = 0.1;
%     parameters(find(strcmp(parameter_names, 'g_BK'))) = 0.02;
%     parameters(find(strcmp(parameter_names, 'g_Kv'))) = 1.0;
%     parameters(find(strcmp(parameter_names, 'g_HERG'))) = 0;
    parameters(find(strcmp(parameter_names, 'V_K'))) = -92;
    
    % Ramp
    v_rest = -80;
    v_pulse = (-100:1:50);
    rest_length = 1000;
    pulse_length = 0.75;
    prepulse_length = 0;
    end_length = rest_length;
    num_pulses = length(v_pulse);
    cycle_length = NaN;
    current = 'ramp';
    plotting = 'true';
    options = odeset('MaxStep',1);

    v_clamp_amps = [v_rest,v_pulse,v_rest];
    v_clamp_times = [rest_length,rest_length+pulse_length:pulse_length:num_pulses*pulse_length+rest_length,num_pulses*pulse_length+rest_length+end_length];
    R_clamp = 0.001; % (GOhm)
    
    init_tspan = [0,10000];
    [T_init, Y_init] = ode15s(rhs_model,init_tspan, states, options, parameters, v_rest, init_tspan(end), R_clamp);
    
    tspan = [0,v_clamp_times(end)];
    [T, Y] = ode15s(rhs_model,tspan, Y_init(end,:), options, parameters, v_clamp_amps, v_clamp_times, R_clamp);

    try
        outs = analyze_v_clamp(Y,T,parameters,state_names, cycle_length, prepulse_length, end_length, rhs_model, v_clamp_amps, v_clamp_times, R_clamp, current, plotting);
        ramp_IK = outs{1};
    catch ME
        reason = {ME.identifier};
        disp(reason)
        ramp_IK = NaN;
    end
end