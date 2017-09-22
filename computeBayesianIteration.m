%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute a Bayesian inference iteration, ie. calculate the posterior probability (and pairwise matching score & likelihood) given the history so far of matching scores, likelihood, priors, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [runningWinScore, runningLikelihood, runningPosterior, assignedMatch] = computeBayesianIteration(runningWinScore, runningLikelihood, runningPosterior, assignedMatch, yCam, yUAV, currT, dims, runningCorrWinSizes, N, M, derivFiltHalfWinSize, sigmaLikelihood, normalizationByRowAndColumn)
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
	if nargin<12 || isempty(derivFiltHalfWinSize)
		derivFiltHalfWinSize = 10;
	end
	if nargin<13 || isempty(sigmaLikelihood)
		sigmaLikelihood = 1.5;
	end
	if nargin<14 || isempty(normalizationByRowAndColumn)
		normalizationByRowAndColumn = 1;
	end

	for iW = 1:length(runningCorrWinSizes)
		winSize = runningCorrWinSizes(iW);
		if mod(currT-1, winSize) == 0  % Check if we just filled out a time-window. If so, compute the iteration
			% Compute the score for this time-window (use a delayed window to avoid using points affected by the derivative filter time-window)
			inds = ((currT-winSize):currT-1)-derivFiltHalfWinSize+1; inds(inds<1) = []; % Indices that represent the points in the current time-window (which lags by derivFilterHalfWinSize points behind currT)
			warning('off', 'stats:pdist2:ConstantPoints'); warning('off', 'stats:pdist2:ZeroPoints');	% Disable correlation and cosine pdist2 warnings
			for iD = 1:length(dims)
				%runningWinScore(:,:,currT,iD,iW) = ones(N,M) - pdist2(yUAV(:,inds,iD)+eps*rand(N,length(inds),1), yCam(:,inds,iD)+eps*rand(M,length(inds),1), 'correlation'); %'cosine');	% pdist2 produces an MxN matrix (it compares each row in X with each row in Y) so we need to transpose it to get NxM
				runningWinScore(:,:,currT,iD,iW) = ones(N,M) - pdist2(yUAV(:,inds,iD)-mean(yUAV(:,inds,iD),2, 'omitnan'), yCam(:,inds,iD)-mean(yCam(:,inds,iD),2, 'omitnan'), 'euclidean')./(repmat(sqrt(sum(yUAV(:,inds,iD).^2,2)), 1,M)+repmat(sqrt(sum(yCam(:,inds,iD).^2,2))', N,1));	% pdist2 produces an MxN matrix (it compares each row in X with each row in Y) so we need to transpose it to get NxM
				%runningWinScore(:,:,currT,iD,iW) = ones(N,M) - pdist2(yUAV(:,inds,iD), yCam(:,inds,iD), 'euclidean')./(repmat(sqrt(sum(yUAV(:,inds,iD).^2,2)), 1,M)+repmat(sqrt(sum(yCam(:,inds,iD).^2,2))', N,1));	% pdist2 produces an MxN matrix (it compares each row in X with each row in Y) so we need to transpose it to get NxM
			end
			
			% Then update the likelihood and posterior, and recompute the optimal id assignment
			[runningLikelihood(:,:,currT,iW), runningPosterior(:,:,currT,iW)] = computeBayesianIterationPosterior(reshape(runningWinScore(:,:,currT,:,iW), N,M,[]), runningPosterior(:,:,currT-winSize,iW), sigmaLikelihood, normalizationByRowAndColumn);
			[assignments, unassignedUAVs, unassignedCams] = assignDetectionsToTracks(-log(runningPosterior(:,:,currT,iW)),-log(0.01));
			assignedMatch(assignments(:,1),currT,iW) = assignments(:,2);
			
			% Finally copy the value of the last iteration up until currT-1
			runningWinScore(:,:,currT-winSize+1:currT-1,:,iW) = repmat(runningWinScore(:,:,currT-winSize,:,iW), 1,1,winSize-1);
			runningLikelihood(:,:,currT-winSize+1:currT-1,iW) = repmat(runningLikelihood(:,:,currT-winSize,iW), 1,1,winSize-1);
			runningPosterior(:,:,currT-winSize+1:currT-1,iW) = repmat(runningPosterior(:,:,currT-winSize,iW), 1,1,winSize-1);
			assignedMatch(:,currT-winSize+1:currT-1,iW) = repmat(assignedMatch(:,currT-winSize,iW), 1,1,winSize-1);
% 		else
% 			% What to do if it's not time to recompute scores yet
% 			runningWinScore(:,:,currT,:,iW) = runningWinScore(:,:,currT-1,:,iW);
% 			runningLikelihood(:,:,currT,iW) = runningLikelihood(:,:,currT-1,iW);
% 			runningPosterior(:,:,currT,iW) = runningPosterior(:,:,currT-1,iW);
% 			assignedMatch(:,currT,iW) = assignedMatch(:,currT-1,iW);
		end
	end
end
