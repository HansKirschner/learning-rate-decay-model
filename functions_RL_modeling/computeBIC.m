function BIC=computeBIC(logLikelihood, numParams, numTrails)

%% computes Bayesian Information Criteria
% When fitting models, it is possible to increase the likelihood by adding parameters, 
% but doing so may result in overfitting. BIC resolves this problem by introducing a penalty term 
% for the number of parameters in the model
% this formuals ist based on Stephan et al., 2009
% where next to the maximum log-likelihood, the number of free parameters, scaled by the
% logarithm of the number of trials is taken into account

% Reference
% Stephan, K. E., Penny, W. D., Daunizeau, J., Moran, R. J., and Friston, K. J. (2009). Bayesian model
% selection for group studies. NeuroImage, 46(4):1004–1017.

BIC=logLikelihood-((numParams/2)*log(numTrails));