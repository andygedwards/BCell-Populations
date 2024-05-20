function modParam_scaling = getScalingFactors(stdev, nModParams, nTrials)
% This function generates random, normally distributed perturbations to
% model parameters and stores the resulting values in the 2D
% 'modParam_scaling' matrix where the rows correspond to trials and the
% columns correspond to the parameters being perturbed.

    baseline = 1;
    sigmaG = stdev.*baseline;
    modParam_scaling = baseline.*exp(sigmaG.*randn(nTrials, nModParams));

end