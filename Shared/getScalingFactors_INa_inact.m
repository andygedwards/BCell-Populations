function modParam_scaling = getScalingFactors_INa_inact(stdev, modParam_names, nTrials)
% This function generates random, log-normally distributed perturbations to
% INa inactivation model parameters and INa activation V_1/2. 
% The resulting parameter values are returned in the 2D
% 'modParam_scaling' matrix where the rows correspond to trials and the
% columns correspond to the parameters being perturbed.
%
% Passed arguments: 
% stdev - a vector of standard deviation values for all parameters to be varied
% nModParams - the number of parameters to be modified, equal to length(stdev)
% modParam_names - the names of the parameters as in the model definition
% nTrials - the number of parameter instances to be generated for each of the nModParams
    n_h_ind = find(contains(modParam_names,'n_h'),1,'first');
    V_m_ind = find(contains(modParam_names,'V_m'),1,'first');
    V_h_ind = find(contains(modParam_names,'V_h'),1,'first');
    mu = 1;
    for i = 1:length(modParam_names)
        sigma = stdev(i);
        norm_mu = log((mu^2) / sqrt(sigma^2 + mu^2));
        norm_sigma = sqrt(log((sigma^2/mu^2) + 1));
        if ~isempty(find(contains(modParam_names(i),'g_Na'),1,'first')) || ~isempty(find(contains(modParam_names(i),'V_h'),1,'first')) % check if the ith parameter name contains g_Na, V_mNa or V_hNa 
            if ~isempty(find(contains(modParam_names(i),'V_h'),1,'first')) % Identify if it is the V_hNa parameter
                samples = lognrnd(norm_mu, norm_sigma, [nTrials, 1]); % Create a log-normal sampling for of length nTrials
                shifted_samples  = mu-(samples-mu); % Ensure that the skewing imposed by log-normal sampling is rightward (toward zero) rather than leftward
                modParam_scaling(:,i) = shifted_samples; % Apply the shift to V_hNa and V_mNa
                modParam_scaling(:,V_m_ind) = shifted_samples;
                modParam_scaling(:,n_h_ind) = mu+(mu-shifted_samples)/1.5; % When shifting V_hNa the steepness of the Hill function must be shifted reciprocally,...
                   % otherwise unphysiologic Hill curves result, and many cells are unstable in V-clamp and I-clamp due to very large window currents
            else
                modParam_scaling(:,i) = lognrnd(norm_mu, norm_sigma, [nTrials, 1]); % if it is g_Na rather than V_hNa, conduct regular log-normal sampling with no shift
            end
        end
        % figure
        % subplot(2,2,1), histogram(modParam_scaling(:,V_h_ind)), title(gca,'V_{hNa}')
        % subplot(2,2,2), histogram(modParam_scaling(:,V_m_ind)), title(gca,'V_{mNa}')
        % subplot(2,2,3), histogram(modParam_scaling(:,n_h_ind)), title(gca,'n_{hNa}')
        % subplot(2,2,4), histogram(modParam_scaling(:,find(contains(modParam_names,'g_Na'),1,'first'))), title(gca,'g_{Na}')
    end
end