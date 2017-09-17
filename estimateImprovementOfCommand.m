%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimates the new posterior probability of an assignment given a command (and predicting how each drone would move based on the assignment)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newP = estimateImprovementOfCommand(command, assignments, runningPrior, defaultScore, sigmaNoiseAccel, sigmaLikelihood)
	if nargin<4 || isempty(defaultScore)
		defaultScore = 0.25;
	end
	if nargin<5 || isempty(sigmaNoiseAccel)
		sigmaNoiseAccel = 0.25;
	end
	if nargin<6
		sigmaLikelihood = [];	% Use default value if not provided
	end
	runningWinScore = defaultScore.*ones(size(runningPrior,1), size(runningPrior,2), 3);
	
	accels = predictPosVelAccelFromCommand(command(3*assignments(:,2)-2), command(3*assignments(:,2)-1), command(3*assignments(:,2)));
	accels = accels + sigmaNoiseAccel.*randn(size(accels));
	for i = 1:size(assignments,1)
		for j = 1:size(assignments,1)
			if j > i
				for iD = 1:3  % Compute the score
					runningWinScore(i,j,iD) = 1-pdist2(accels(iD,:,i), accels(iD,:,j), 'euclidean')/(eps+norm(accels(iD,:,i))+norm(accels(iD,:,j)));
				end
			elseif j < i  % Already computed
				runningWinScore(i,j,:) = runningWinScore(j,i,:);
			else % j == i  % You with yourself should have a high score
				runningWinScore(i,i,:) = 0.9;
			end
		end
	end
	
	% Based on these estimated scores, compute a bayesian iteration and compute the new probability of such assignment
	[~, runningPosterior] = computeBayesianIterationPosterior(runningWinScore, runningPrior, sigmaLikelihood);
	newP = prod(runningPosterior(sub2ind(size(runningPosterior), assignments(:,1),assignments(:,2))));
end
