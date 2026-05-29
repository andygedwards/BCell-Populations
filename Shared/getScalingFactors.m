function modParam_scaling = getScalingFactors(stdev, nModParams, nTrials, varargin)
% This function generates random, normally distributed perturbations to
% model parameters and stores the resulting values in the 2D
% 'modParam_scaling' matrix where the rows correspond to trials and the
% columns correspond to the parameters being perturbed.

    mu = 1;
    for i = 1:nModParams
        sigma = stdev(i);
        norm_mu = log((mu^2) / sqrt(sigma^2 + mu^2));
        norm_sigma = sqrt(log((sigma^2/mu^2) + 1));
        modParam_scaling(:,i) = lognrnd(norm_mu, norm_sigma, [nTrials, 1]);
    end
    if ~isempty(varargin)
        PQcol = varargin{1};
        thresh = 8.25;
        initHighInds = find(modParam_scaling(:,PQcol)>thresh);
        while ~isempty(find(modParam_scaling(:,PQcol)>thresh))
            highInds = find(modParam_scaling(:,PQcol)>thresh);
            corrected = abs(thresh-(modParam_scaling(highInds,PQcol)-thresh));
            modParam_scaling(highInds,PQcol) = corrected;
        end
    end
end