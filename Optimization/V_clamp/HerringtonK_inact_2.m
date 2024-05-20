function [Amp1_IK_inact, tau1_IK_inact, Amp2_IK_inact, tau2_IK_inact] = HerringtonK_inact_2(parameters, parameter_names, states, state_names, rhs_model)
    % Herrington protocols
%     parameters(find(strcmp(parameter_names, 'g_KATP0'))) = 0.01;
%     parameters(find(strcmp(parameter_names, 'g_KATP_hat'))) = 0.05;
%     parameters(find(strcmp(parameter_names, 'g_SK'))) = 0.1;
%     parameters(find(strcmp(parameter_names, 'g_BK'))) = 0.02;
%     parameters(find(strcmp(parameter_names, 'g_Kv'))) = 1.0;
%     parameters(find(strcmp(parameter_names, 'g_HERG'))) = 0;
    parameters(find(strcmp(parameter_names, 'V_K'))) = -92;
    
    % K_inact_2

    v_rest = -80;
    v_pulse = 20;
    rest_length = 10000;
    pulse_length = 10000;
    cycle_length = rest_length+pulse_length;
    prepulse_length = 0;
    end_length = 0;
    num_pulses = length(v_rest);
    current = 'K_inact_2';
    plotting = 'true';
    options = odeset('MaxStep',1);

    v_clamp_amps(1:2:2*num_pulses) = v_rest; v_clamp_amps(2:2:num_pulses*2) = v_pulse; 
    v_clamp_times(1:2:2*num_pulses) = rest_length:cycle_length:num_pulses*cycle_length; v_clamp_times(2:2:num_pulses*2) = cycle_length:cycle_length:num_pulses*cycle_length; 
    R_clamp = 0.001; % (GOhm)
  
    init_tspan = [0,10000];
    [T_init, Y_init] = ode15s(rhs_model,init_tspan, states, options, parameters, v_rest, init_tspan(end), R_clamp);
    
    tspan = [0,v_clamp_times(end)];
    [T, Y] = ode15s(rhs_model,tspan, Y_init(end,:), options, parameters, v_clamp_amps, v_clamp_times, R_clamp);

    try
        outs = analyze_v_clamp(Y,T,parameters,state_names, cycle_length, prepulse_length, end_length, rhs_model, v_clamp_amps, v_clamp_times, R_clamp, current, plotting);
        Amp1_IK_inact = outs{1};
        tau1_IK_inact = outs{2};
        Amp2_IK_inact = outs{3};
        tau2_IK_inact = outs{4};
    catch ME
        reason = {ME.identifier};
        disp(reason)
        Amp1_IK_inact = NaN;
        tau1_IK_inact = NaN;
        Amp2_IK_inact = NaN;
        tau2_IK_inact = NaN;
    end
end