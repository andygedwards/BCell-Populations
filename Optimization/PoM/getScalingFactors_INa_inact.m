function modParam_scaling = getScalingFactors_INa_inact(stdev, modParam_names, nTrials)
% This function generates random, normally distributed perturbations to
% INa inactivation model parameters and INa activation V_1/2. The variation in kinetic parameters is
% coordinated such that activation and inactivation curves are altered in an approximately corresponding manner 
% (i.e. right-shifted or left-shifted to similar degrees). The resulting parameter values are returned in the 2D
% 'modParam_scaling' matrix where the rows correspond to trials and the
% columns correspond to the parameters being perturbed.
%
% Passed arguments: 
% stdev - a vector of standard deviation values for all parameters to be varied
% nModParams - the number of parameters to be modified, equal to length(stdev)
% modParam_names - the names of the parameters as in the model definition
% nTrials - the number of parameter instances to be generated for each of the nModParams

    baseline = 1;
    sigmaG = stdev.*baseline;
    gNa_ind = find(contains(modParam_names,'g_Na'),1,'first');
    V_h_ind = find(contains(modParam_names,'V_h'),1,'first');
    n_h_ind = find(contains(modParam_names,'n_h'),1,'first');
    V_m_ind = find(contains(modParam_names,'V_m'),1,'first');
    modParam_scaling(:,gNa_ind) = baseline*exp(sigmaG(gNa_ind)*randn(nTrials, length(gNa_ind)));
    modParam_scaling(:,V_h_ind) = baseline*exp(sigmaG(V_h_ind)*randn(nTrials, length(V_h_ind)));
    modParam_scaling(:,V_m_ind) = baseline+(baseline-modParam_scaling(:,V_h_ind))*1.5;
    modParam_scaling(:,n_h_ind) = baseline+(baseline-modParam_scaling(:,V_h_ind))/1.5;
end