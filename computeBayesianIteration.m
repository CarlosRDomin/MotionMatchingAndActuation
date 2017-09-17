%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute a Bayesian inference iteration, ie. calculate the posterior probability (and pairwise matching score & likelihood) given the history so far of matching scores, likelihood, priors, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [runningWinScore, runningLikelihood, runningPrior, assignedMatch] = computeBayesianIteration(runningWinScore, runningLikelihood, runningPrior, assignedMatch, yCam, yUAV, currT, dims, runningCorrWinSizes, N, M, sigmaLikelihood)
	if nargin<8 || isempty(dims)
		dims = 1:3;
	end
	if nargin<9 || isempty(runningCorrWinSizes)
		runningCorrWinSizes = 15;
	end
	if nargin<10 || isempty(N)
		N = size(yUAV,1);
	end
	if nargin<11 || isempty(M)
		M = size(yCam,1);
	end
	if nargin<12 || isempty(sigmaLikelihood)
		sigmaLikelihood = 1.5;
	end

	for iW = 1:length(runningCorrWinSizes)
		winSize = runningCorrWinSizes(iW);
		if mod(currT, winSize) == 0  % Check if we just filled out a time-window. If so, compute the iteration
			inds = (currT-winSize+1):currT; % Indices that represent the points in the current time-window
			for iD = 1:length(dims)
				warning('off', 'stats:pdist2:ConstantPoints'); warning('off', 'stats:pdist2:ZeroPoints');	% Disable correlation and cosine pdist2 warnings
				%runningWinScore(:,:,currT,iD,iW) = ones(N,M) - pdist2(yUAV(:,inds,iD)+eps*rand(N,length(inds),1), yCam(:,inds,iD)+eps*rand(M,length(inds),1), 'correlation'); %'cosine');	% pdist2 produces an MxN matrix (it compares each row in X with each row in Y) so we need to transpose it to get NxM
				runningWinScore(:,:,currT,iD,iW) = ones(N,M) - pdist2(yUAV(:,inds,iD)-mean(yUAV(:,inds,iD),2, 'omitnan'), yCam(:,inds,iD)-mean(yCam(:,inds,iD),2, 'omitnan'), 'euclidean')./(repmat(sqrt(sum(yUAV(:,inds,iD).^2,2)), 1,M)+repmat(sqrt(sum(yCam(:,inds,iD).^2,2))', N,1));	% pdist2 produces an MxN matrix (it compares each row in X with each row in Y) so we need to transpose it to get NxM
				%runningWinScore(:,:,currT,iD,iW) = ones(N,M) - pdist2(yUAV(:,inds,iD), yCam(:,inds,iD), 'euclidean')./(repmat(sqrt(sum(yUAV(:,inds,iD).^2,2)), 1,M)+repmat(sqrt(sum(yCam(:,inds,iD).^2,2))', N,1));	% pdist2 produces an MxN matrix (it compares each row in X with each row in Y) so we need to transpose it to get NxM
			end
			
			[runningLikelihood(:,:,currT,iW), runningPrior(:,:,currT+1,iW)] = computeBayesianIterationPosterior(reshape(runningWinScore(:,:,currT,:,iW), N,M,[]), runningPrior(:,:,currT+1-winSize,iW), sigmaLikelihood);
			[assignments, unassignedUAVs, unassignedCams] = assignDetectionsToTracks(-log(runningPrior(:,:,currT+1,iW)),-log(0.01));
			assignedMatch(assignments(:,1),currT,iW) = assignments(:,2);
		else
			% What to do if it's not time to recompute scores yet
			for iD = 1:length(dims)
				runningWinScore(:,:,currT,iD,iW) = runningWinScore(:,:,currT-1,iD,iW);
			end
			runningLikelihood(:,:,currT,iW) = runningLikelihood(:,:,currT-1,iW);
			runningPrior(:,:,currT+1,iW) = runningPrior(:,:,currT,iW);
			assignedMatch(:,currT,iW) = assignedMatch(:,currT-1,iW);
		end
	end
end
