close all; clear all;
cdToWorkingDirectory;
includeScripts;

N = 5;
M = N;
dims = 1;
pointsPerIter = 100;
motionAccel = zeros(3, pointsPerIter, length(dims));
motionVel = zeros(3, pointsPerIter, length(dims));
motionPos = zeros(3, pointsPerIter, length(dims));
deltaP = 1;	% 1m per timestep
deltaT = 1;	% Timestep/iteration: 1s
K = 2*pi*deltaP/deltaT.^2; % a(t) = K*sin(2*pi*(1/deltaT)*t)
t = deltaT/pointsPerIter*(1:pointsPerIter); % t = 0:deltaT/pointsPerIter:deltaT;
tt = 0:deltaT/(10*pointsPerIter):deltaT; %tt=deltaT*(0:0.5:1);
for i = -1:1
	aa = K*i*sin(2*pi/deltaT*tt); aa(end)=0; % Improve precision
	vv = cumtrapz(tt,aa); vv(end)=0; % Improve precision
	pp = cumtrapz(tt,vv);
	a = interp1(tt,aa, t);
	v = interp1(tt,vv, t);
	p = interp1(tt,pp, t);
	
% 	vv=[0, 2*i*deltaP/deltaT, 0];
% 	v=interp1(tt,vv, t);
% 	p=cumtrapz(t,v);
% % 	a=diff([v(1) v(1:end-1); v(2:end) v(end)])./diff([t(1) t(1:end-1); t(2:end) t(end)]);
% % 	a=diff([v v(end); v(1) v])./diff([t t(end); t(1) t]);
% 	a=[diff(v)./diff(t) 0];
	motionAccel(i+2,:,:) = a;
	motionVel(i+2,:,:) = v;
	motionPos(i+2,:,:) = p;
end
%%
tMax = 10; t = 0:deltaT/pointsPerIter:tMax; runningCorrWinSizes = pointsPerIter;
yCam = cat(2, zeros(M, 1, length(dims)), NaN(M, length(t)-1, length(dims)));
yUAV = cat(2, zeros(N, 1, length(dims)), NaN(N, length(t)-1, length(dims)));
posUAV = cat(2, 10*randn(N, 1, length(dims)), NaN(N, length(t)-1, length(dims)));
runningCorr = cat(3, zeros(N, M, 1, length(dims), length(runningCorrWinSizes)), NaN(N, M, length(t)-1, length(dims), length(runningCorrWinSizes)));
assignedMatch = NaN(M, length(t), length(dims), length(runningCorrWinSizes));
muSgivenTheta = eye(N); sigmaSgivenTheta = repmat(1.5*eye(N), 1,1,N);
runningPrior = cat(3, ones(N, M, 1, length(dims), length(runningCorrWinSizes))./N, NaN(N, M, length(t), length(dims), length(runningCorrWinSizes)));
runningLikelihood = NaN(N, M, length(t), length(dims), length(runningCorrWinSizes));
groundTruthAssignment = randperm(N,M); groundTruthAssignment=1:M;	% Get a random permutation of M elements picked from the set 1:N
sigmaNoise = 1;	% m/s2

figure('Units','normalized', 'Position',[0.3 0.4 0.4 0.2]); colormap('jet'); scatter(posUAV(:,1,:), 0*posUAV(:,1,:), [],1:N, 'filled'); xlim([-25, 25]); title('Starting pos');
%%
for currT = 0:deltaT:(tMax-deltaT)
	currTind = currT*pointsPerIter/deltaT + 1;
	risk = zeros(1, 3^N);
	uniqueness = zeros(1, 3^N);
	for i = 1:3^N
		motionMagns = permn(-1:1, N,i); %decomposeMotion(i,N).*deltaP;
		confidence = prod(runningPrior(:,:,currTind,:,1), 4);	% Confidence = P_x * P_y * P_z -> Multiply at dimension 4 (dims)
		yUAV(:,currTind+(1:pointsPerIter),:) = motionAccel(motionMagns+2,:,:);
		posUAV(:,currTind+(1:pointsPerIter),:) = posUAV(:,currTind,:) + motionPos(motionMagns+2,:,:);
% 		expectedPos = pos;
% 		for j = 1:n
% 			expectedPos(j) = expectedPos(j) + motionMagns*confidence(:,j);
% 		end

		P_perm_total = 0;
		combinations = nchoosek(1:N, M);				% Get all the possible subgroups of M (out of N) drones
		for j = 1:size(combinations,1)
			permutations = perms(combinations(j,:));	% Get all the possible permutations within the selected combination (eg: all permutations of [1 2 4] being assigned to objects 1:M)
			for k = 1:size(permutations,1)
				P_perm = prod(confidence(sub2ind(size(confidence), permutations(k,:), 1:M))); P_perm_total = P_perm_total + P_perm;
				yCam(:,currTind+(1:pointsPerIter),:) = motionAccel(motionMagns(permutations(k,:))+2,:,:);
				[runningCorr, runningLikelihood, runningPrior, assignedMatch] = computeBayesianIteration(runningCorr, runningLikelihood, runningPrior, assignedMatch, yCam, yUAV, currTind+pointsPerIter-1, dims, runningCorrWinSizes, N, M, muSgivenTheta, sigmaSgivenTheta);
				risk(i) = risk(i) + P_perm*min(pdist(posUAV(:,currTind+pointsPerIter,:)));
				new_confidence = prod(runningPrior(:,:,currTind+pointsPerIter,:,1), 4);
				%uniqueness(i) = uniqueness(i) + P_perm*prod(new_confidence(sub2ind(size(new_confidence), permutations(k,:), 1:M)))/prod(prod(new_confidence));
				%uniqueness(i) = uniqueness(i) + P_perm*sum(new_confidence(sub2ind(size(new_confidence), permutations(k,:), 1:M)));
				%aux = new_confidence(:)./sum(new_confidence(:)); aux(aux==0)=1; % Avoid errors taking 0.*log2(0) (which should give 0, like 1.*log2(1))
				%uniqueness(i) = uniqueness(i) + P_perm*sum(aux.*log2(aux));
				uniqueness(i) = uniqueness(i) - P_perm*1e6*prod(new_confidence(:));
				dispImproved(sprintf('@t=%2d - Trying permutation %3d/%3d for motion command %5d/%5d\t(%6.2f%%)', currT, k,size(permutations,1), i,3^N, 100*((i-1)*size(permutations,1)+k)/(size(permutations,1)*3^N)));
			end
		end
		risk(i) = risk(i)/P_perm_total;
		uniqueness(i) = uniqueness(i)/P_perm_total;
% 		expectedPos = posUAV(:,currT,:) + (motionMagns*runningPrior(:,1:N))';
% 		risk(i) = min(pdist(expectedPos));
	end
	
	[riskVal, riskInd] = max(risk);
	[uniquenessVal, uniquenessInd] = max(uniqueness);
	motionMagns = permn(-1:1, N,uniquenessInd);
	dispImproved(sprintf('@t=%2d - Sending motion command %s\n', currT, mat2str(motionMagns)), 'keepthis');
	yCam(:,currTind+(1:pointsPerIter),:) = motionAccel(motionMagns(groundTruthAssignment)+2,:,:)+sigmaNoise.*randn(M,pointsPerIter,length(dims));
	yUAV(:,currTind+(1:pointsPerIter),:) = motionAccel(motionMagns+2,:,:)+sigmaNoise.*randn(N,pointsPerIter,length(dims));
	posUAV(:,currTind+(1:pointsPerIter),:) = posUAV(:,currTind,:) + motionPos(motionMagns+2,:,:);
	[runningCorr, runningLikelihood, runningPrior, assignedMatch] = computeBayesianIteration(runningCorr, runningLikelihood, runningPrior, assignedMatch, yCam, yUAV, currTind+pointsPerIter-1, dims, runningCorrWinSizes, N, M, muSgivenTheta, sigmaSgivenTheta);
	figure('Units','normalized', 'Position',[0.3 0.4 0.4 0.2]); colormap('jet'); hold on;
	scatter(posUAV(:,currTind+pointsPerIter,:), 0*posUAV(:,currTind+pointsPerIter,:), [],1:N, 'filled'); xlim([-25, 25]); title(sprintf('At t=%d, motion:%s', currT+deltaT, mat2str(motionMagns)));
end
%%
runningPrior(:,:,2,:,:) = runningPrior(:,:,1,:,:);
for currT = 2:length(t)-1	% Fill in the gaps
	[runningCorr, runningLikelihood, runningPrior, assignedMatch] = computeBayesianIteration(runningCorr, runningLikelihood, runningPrior, assignedMatch, yCam, yUAV, currT, dims, runningCorrWinSizes, N, M, muSgivenTheta, sigmaSgivenTheta);
end
%%
outFields = {'runningCorr','runningLikelihood','runningPrior','assignedMatch','N','M','iCams','dims','runningCorrWinSizes','cropToMinT','t','yCam','yUAV'};
runningCorrStruct = cell2struct(cell(1, length(outFields)), outFields, 2);
iCams=1:M; cropToMinT=false;
for f = outFields	% Populate output struct with results
	runningCorrStruct.(f{:}) = eval(f{:});
end
plotExperimentResults(runningCorrStruct, struct('rawAccel',true, 'scatterCorr',false, 'runningLikelihoodVsWinSize',true, 'runningLikelihoodFull',true));