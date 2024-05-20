function [early_exocytosis, late_exocytosis, total_exocytosis, v_half, vh_r2, peak_INa, v_half_act, vha_r2, peak_ICa, late_ICa]... , v_half_IK_inact], tau1_IK_inact, tau2_IK_inact, Amp1_IK_inact, Amp2_IK_inact]...
    = run_v_clamp_sens(parameters, parameter_names, states, state_names, rhs_model)

% Run all the voltage clamp modules

    % @ANDY --> check the HerringtonRamp and HerringtonK_square protocols >
    % weird behaviour.
                      
    [v_half, vh_r2] = INa_inact(parameters, parameter_names, states, state_names, rhs_model);
    [peak_INa, v_half_act, vha_r2, peak_ICa, late_ICa] = PeakCurrents(parameters, parameter_names, states, state_names, rhs_model);
    [early_exocytosis, late_exocytosis, total_exocytosis] = Exocytosis(parameters, parameter_names, states, state_names, rhs_model);
    %[ramp_IK] = HerringtonRamp(parameters, parameter_names, states, state_names, rhs_model);
    %[Peak_IK_minus20, Peak_IK_max] = HerringtonK_square_sens(parameters, parameter_names, states, state_names, rhs_model);
    %[v_half_act_IK, n_act_IK] = HerringtonK_tail_1_sens(parameters, parameter_names, states, state_names, rhs_model);
    %[v_half_IK_inact] = HerringtonK_inact_1(parameters, parameter_names, states, state_names, rhs_model)
    %[Amp1_IK_inact, tau1_IK_inact, Amp2_IK_inact, tau2_IK_inact] = HerringtonK_inact_2(parameters, parameter_names, states, state_names, rhs_model)

end