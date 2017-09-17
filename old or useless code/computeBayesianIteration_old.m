function [runningCorr, runningLikelihood, runningPrior, assignedMatch] = computeBayesianIteration(runningCorr, runningLikelihood, runningPrior, assignedMatch, yCam, yUAV, currT, dims, runningCorrWinSizes, N, M, muSgivenTheta, sigmaSgivenTheta)
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
	if nargin<12 || isempty(muSgivenTheta)
		muSgivenTheta = eye(N);
	end
	if nargin<13 || isempty(sigmaSgivenTheta)
		sigmaSgivenTheta = repmat(1.5*eye(N), 1,1,N);
	end

	for iW = 1:length(runningCorrWinSizes)
		winSize = runningCorrWinSizes(iW);
		if mod(currT, winSize) == 0
			inds = (currT-winSize+1):currT;
			for iD = 1:length(dims)
				warning('off', 'stats:pdist2:ConstantPoints'); warning('off', 'stats:pdist2:ZeroPoints');	% Disable correlation and cosine pdist2 warnings
				runningCorr(:,:,currT,iD,iW) = ones(N,M) - pdist2(yUAV(:,inds,iD)+eps*rand(size(yCam,1),length(inds),1), yCam(:,inds,iD)+eps*rand(size(yCam,1),length(inds),1), 'correlation'); %'cosine');	% pdist2 produces an MxN matrix (it compares each row in X with each row in Y) so we need to transpose it to get NxM
				%runningCorr(:,:,currT,iD,iW) = ones(N,M) - pdist2(yUAV(:,inds,iD)-mean(yUAV(:,inds,iD),2), yCam(:,inds,iD)-mean(yCam(:,inds,iD),2), 'euclidean')./(repmat(sqrt(sum(yUAV(:,inds,iD).^2,2)), 1,M)+repmat(sqrt(sum(yCam(:,inds,iD).^2,2))', N,1));	% pdist2 produces an MxN matrix (it compares each row in X with each row in Y) so we need to transpose it to get NxM
				%runningCorr(:,:,currT,iD,iW) = ones(N,M) - pdist2(yUAV(:,inds,iD), yCam(:,inds,iD), 'euclidean')./(repmat(sqrt(sum(yUAV(:,inds,iD).^2,2)), 1,M)+repmat(sqrt(sum(yCam(:,inds,iD).^2,2))', N,1));	% pdist2 produces an MxN matrix (it compares each row in X with each row in Y) so we need to transpose it to get NxM
				
				% From the correlation, update the likelihood of each column (object iC in the cam vs every UAV). 
% 				for iC = find(sum(isnan(runningCorr(:,:,currT,iD,iW)))<N)
% 					c = runningCorr(:,iC,currT,iD,iW)';
% 					validUAVs = find(~isnan(c));	% Make sure to avoid NaN values (which happen when a yUAV is shorter in time than the yCam, so pdist2 can't evaluate the distance)
% 					runningLikelihood(validUAVs,iC,currT,iD,iW) = mvnpdf(c(validUAVs), muSgivenTheta(validUAVs,validUAVs), sigmaSgivenTheta(validUAVs,validUAVs,validUAVs)); % runningLikelihood(:,iCam,endWin,d,iW) = runningLikelihood(:,iCam,endWin,d,iW)./sum(runningLikelihood(:,iCam,endWin,d,iW));
% 				end
				runningLikelihood(:,:,currT,iD,iW) = normpdf(runningCorr(:,:,currT,iD,iW), 1,2/3);
				%runningLikelihood(:,:,currT,iD,iW) = runningCorr(:,:,currT,iD,iW);
				
				runningPrior(:,:,currT+1,iD,iW) = runningLikelihood(:,:,currT,iD,iW).*runningPrior(:,:,currT+1-winSize,iD,iW);
				runningPrior(:,:,currT+1,iD,iW) = runningPrior(:,:,currT+1,iD,iW)./(repmat(sum(runningPrior(:,:,currT+1,iD,iW),1,'omitNaN'), N,1)+repmat(sum(runningPrior(:,:,currT+1,iD,iW),2,'omitNaN'), 1,M)-runningPrior(:,:,currT+1,iD,iW));
				
				% for iC = 1:M
				% 	if isequal(runningPrior(:,iC,endWin+1,iD,iW), zeros(N,1)), runningPrior(:,iC,endWin+1,iD,iW) = 1/N; end
				% end
% 				[m, assignedMatch(:,currT,iD,iW)] = max(runningPrior(:,:,currT+1,iD,iW));
% 				assignedMatch(isnan(m),currT,iD,iW) = NaN;
				
				[assignments, unassignedUAVs, unassignedCams] = assignDetectionsToTracks(-log(runningPrior(:,:,currT+1,iD,iW)),-log(0.01));
				assignedMatch(assignments(:,1),currT,iD,iW) = assignments(:,2);
			end
		else
			% What to do if it's not time to recompute scores yet
			for iD = 1:length(dims)
				runningCorr(:,:,currT,iD,iW) = runningCorr(:,:,currT-1,iD,iW);
				runningLikelihood(:,:,currT,iD,iW) = runningLikelihood(:,:,currT-1,iD,iW);
				runningPrior(:,:,currT+1,iD,iW) = runningPrior(:,:,currT,iD,iW);
				assignedMatch(:,currT,iD,iW) = assignedMatch(:,currT-1,iD,iW);
			end
		end
	end
end
