%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializes the time-varying variables of our matching framework: accelUAV, accelCam, runningWinScore, runningLikelihood, runningPosterior, assignedMatch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [accelUAV, accelCam, runningWinScore, runningLikelihood, runningPosterior, assignedMatch] = initFrameworkVars(N, M, dims, t, frameworkWinSize, derivFiltHalfWinSize, normalizationByRowAndColumn)
	if nargin<6 || isempty(derivFiltHalfWinSize), derivFiltHalfWinSize = 21; end
	if nargin<7 || isempty(normalizationByRowAndColumn), normalizationByRowAndColumn = 0; end % Regular Bayesian normalization (divide each cell by the sum of its row)
	
	accelUAV = cat(2, zeros(N, 1, length(dims)), NaN(N, length(t)-1, length(dims)));
	accelCam = cat(2, zeros(M, derivFiltHalfWinSize, length(dims)), NaN(M, length(t)-derivFiltHalfWinSize, length(dims)));
	runningWinScore = cat(3, zeros(N, M, 1, length(dims), length(frameworkWinSize)), NaN(N, M, length(t)-1, length(dims), length(frameworkWinSize)));
	runningLikelihood = NaN(N, M, length(t), length(frameworkWinSize));
	if normalizationByRowAndColumn==1, normalizeBy = N+M+1; else, normalizeBy = N; end
	runningPosterior = cat(3, ones(N, M, 1, length(frameworkWinSize))./normalizeBy, NaN(N, M, length(t)-1, length(frameworkWinSize)));
	assignedMatch = cat(2, repmat(randperm(N)', 1,1,length(frameworkWinSize)), NaN(N, length(t)-1, length(frameworkWinSize)));
end

