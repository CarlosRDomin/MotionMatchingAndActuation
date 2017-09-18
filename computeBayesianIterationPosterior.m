%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Bayesian posterior (one of the steps of an inference iteration), ie. calculate the posterior probability (and pairwise matching likelihood) given this iteration's window score and prior probability distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [runningLikelihood, runningPosterior] = computeBayesianIterationPosterior(runningWinScore, runningPrior, sigmaLikelihood)
	%%% This function is part of the computeBayesianIteration process, but I set it up as a function itself so we can only call this part when trying to find the optimal command to send in the next actuation iteration
	if nargin<3 || isempty(sigmaLikelihood)
		sigmaLikelihood = 1;
	end
	N=size(runningWinScore,1);
	M=size(runningWinScore,2);
	
	runningLikelihood = reshape(mvnpdf(reshape(runningWinScore, N*M,[]), ones(1,size(runningWinScore,3)), sigmaLikelihood.*eye(size(runningWinScore,3))), N,M);
	runningLikelihood = runningLikelihood./(repmat(sum(runningLikelihood,1, 'omitNaN'), N,1)+repmat(sum(runningLikelihood,2, 'omitNaN'), 1,M)-runningLikelihood);
	%runningLikelihood = prod(runningWinScore,3);	

	runningPosterior = runningLikelihood.*runningPrior;
	%runningPosterior = runningPosterior./(repmat(sum(runningPosterior,1, 'omitNaN'), N,1)+repmat(sum(runningPosterior,2, 'omitNaN'), 1,M)-runningPosterior);
	runningPosterior = runningPosterior ./ repmat(sum(runningPosterior,2, 'omitNaN'), 1,M);
	runningPosterior(isnan(runningPosterior)) = 0;
	% for iC = 1:M
	% 	if isequal(runningPosterior(:,iC), zeros(N,1)), runningPosterior(:,iC) = 1/N; end
	% end
end
