%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines whether a simulation experiment should end: checks if all N highest posterior values are above a given threshold (eg, 99.99%)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shouldEnd = testIfExperimentShouldEnd(runningPrior, threshAllPosteriors)
	% Get the highest posterior for each column and sort them from higher to lower
	m = sort(max(squeeze(runningPrior), [], 1), 'descend');
	% Then, make sure that the first M (or N, whichever is smaller) values are all above threshAllPosteriors
	minDim = min(size(runningPrior,1), size(runningPrior,2));
	shouldEnd = (m(minDim) >= threshAllPosteriors);
end

