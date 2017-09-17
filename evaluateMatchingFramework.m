%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs a set of experiments with our matching framework and tries to determine which drone is which
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;
cdToThisScriptsDirectory;

dims = 1:3;		% 3-axis acceleration
pointsPerIter = 30;
maxSpeed = 1;	% m/s max speed drones can fly at
deltaT = 1;		% Timestep/iteration: 1s
deltaP = maxSpeed*deltaT;	% Max distance drones can move during a timestep
sigmaNoiseCam = 0.05;	% m, noise in camera's position estimation
sigmaNoisePos = 0.015;	% m, error in position after performing a motion command
sigmaNoiseAccel = 0.25;	% m/s2, noise in IMU's accelerometer data
roomDimensions = [5, 5, 2.5];  % Width x Depth x Height of the room in m
spotterCam = struct('FOV',120*pi/180, 'pos',[], 'orient',[1 0 0; 0 0 -1; 0 1 0]);	% Camera's Field Of View in rad (120º), orientation pointing in the direction of the y-axis
spotterCam.pos = [roomDimensions(1)/2, -(roomDimensions(1)/2) / tan(spotterCam.FOV/2), roomDimensions(3)/2];

scatterPtSize = 80;
scatterCoordLim = 5;
scatterView = [-72, 14];	% degrees of azimuth (theta) and elevation (pi/2-phi)

for N = [5]
	M = N;
	for typeOfMotion = {'landed'} %{'random', 'hovering', 'landed'}
		tMax = 10; t = 0:deltaT/pointsPerIter:tMax; runningCorrWinSizes = pointsPerIter; iW=1;
		accelUAV = cat(2, zeros(N, 1, length(dims)), NaN(N, length(t)-1, length(dims)));
		accelCam = cat(2, zeros(M, 1, length(dims)), NaN(M, length(t)-1, length(dims)));
		if strcmp(typeOfMotion, 'landed')
			initPosUAV = cat(3, rand(N,1, length(dims)-1) .* reshape(roomDimensions(1:length(dims)-1), 1,1,[]), zeros(N,1,1));
		elseif strcmp(typeOfMotion, 'hovering')
			initPosUAV = cat(3, rand(N,1, length(dims)-1) .* reshape(roomDimensions(1:length(dims)-1), 1,1,[]), 0.5*roomDimensions(end)*ones(N,1,1));
		elseif strcmp(typeOfMotion, 'random')
			initPosUAV = rand(N,1, length(dims)) .* reshape(roomDimensions, 1,1,[]);
		end
		posUAVgt = cat(2, initPosUAV, NaN(N, length(t)-1, length(dims)));
		posUAVcam = posUAVgt;
		runningWinScore = cat(3, zeros(N, M, 1, length(dims), length(runningCorrWinSizes)), NaN(N, M, length(t)-1, length(dims), length(runningCorrWinSizes)));
		runningLikelihood = NaN(N, M, length(t), length(runningCorrWinSizes));
		runningPrior = cat(3, ones(N, M, 1, length(runningCorrWinSizes))./(N+M-1), NaN(N, M, length(t), length(runningCorrWinSizes)));
		assignedMatch = NaN(N, length(t), length(runningCorrWinSizes));
		groundTruthAssignment = randperm(N,M)'; groundTruthAssignment=(1:M)';	% Get a random permutation of M elements picked from the set 1:N

		figure('Units','normalized', 'Position',[0.3 0.4 0.4 0.25]);
		plotDronesInRoom(posUAVgt(:,1,:), roomDimensions, spotterCam);
		
		for currT = deltaT*(0:tMax-1)
			currTind = currT*pointsPerIter/deltaT + 1;
			posUAVcam(:,currTind,:) = posUAVgt(:,currTind,:) + sigmaNoiseCam.*randn(size(posUAVcam(:,currTind,:)));

			if strcmp(typeOfMotion, 'landed')
				a = zeros(N, pointsPerIter, length(dims));
				p = zeros(N, pointsPerIter, length(dims));
				posNoise = 0;
			elseif strcmp(typeOfMotion, 'hovering')
				posNoise = sigmaNoisePos;
				command = zeros(3*N,1);
				[a,~,p] = predictPosVelAccelFromCommand(command(1:3:end), command(2:3:end), command(3:3:end)); a=permute(a,3:-1:1); p=permute(p,3:-1:1);
			elseif strcmp(typeOfMotion, 'random')
				posNoise = sigmaNoisePos;
				command = reshape([deltaP*rand(1,N); [2*pi; pi].*rand(2,N)], [],1);
				[a,~,p] = predictPosVelAccelFromCommand(command(1:3:end), command(2:3:end), command(3:3:end)); a=permute(a,3:-1:1); p=permute(p,3:-1:1);
			end
			posUAVgt(:,currTind+(1:pointsPerIter),:) = posUAVgt(:,currTind,:) + posNoise.*randn(N,pointsPerIter,length(dims));
			accelUAV(:,currTind+(1:pointsPerIter),:) = sigmaNoiseAccel.*randn(N,pointsPerIter,length(dims));
			accelCam(:,currTind+(1:pointsPerIter),:) = sigmaNoiseAccel.*randn(M,pointsPerIter,length(dims));
			posUAVgt(groundTruthAssignment,currTind+(1:pointsPerIter),:) = posUAVgt(groundTruthAssignment,currTind+(1:pointsPerIter),:) +p(groundTruthAssignment,:,:);
			accelUAV(:,currTind+(1:pointsPerIter),:) = accelUAV(:,currTind+(1:pointsPerIter),:) + a;
			accelCam(:,currTind+(1:pointsPerIter),:) = accelCam(:,currTind+(1:pointsPerIter),:) + a(groundTruthAssignment,:,:);
			[runningWinScore, runningLikelihood, runningPrior, assignedMatch] = computeBayesianIteration(runningWinScore, runningLikelihood, runningPrior, assignedMatch, accelCam, accelUAV, currTind+pointsPerIter-1, dims, runningCorrWinSizes, N, M);
			dispImproved(sprintf('\nPosterior likelihood:\n%s%s\n', [num2str(100.*runningPrior(:,:,currTind+pointsPerIter,iW), '%8.2f')'; repmat(13,1,N)], repmat('-',1,50)), 'keepthis');

			plotDronesInRoom(posUAVgt(:,currTind+pointsPerIter,:), roomDimensions, spotterCam);
		end
		
		%% Fill in the gaps (simulation only computes everything every pointsPerIter points)
		runningPrior(:,:,2,:) = runningPrior(:,:,1,:);
		for currT = 2:length(t)-1	% Fill in the gaps
			[runningWinScore, runningLikelihood, runningPrior, assignedMatch] = computeBayesianIteration(runningWinScore, runningLikelihood, runningPrior, assignedMatch, accelCam, accelUAV, currT, dims, runningCorrWinSizes, N, M);
		end
		
		%% Plot the experiment results (raw accels, likelihood, posteriors...)
		outFields = {'runningWinScore','runningLikelihood','runningPrior','assignedMatch','N','M','iCams','dims','runningCorrWinSizes','cropToMinT','t','accelCam','accelUAV','posUAVgt','groundTruthAssignment','pointsPerIter','deltaP','deltaT','sigmaNoiseCam','sigmaNoisePos','sigmaNoiseAccel'};
		runningCorrStruct = cell2struct(cell(1, length(outFields)), outFields, 2);
		iCams=1:M; cropToMinT=false;
		for f = outFields	% Populate output struct with results
			runningCorrStruct.(f{:}) = eval(f{:});
		end
		plotExperimentResults(runningCorrStruct, struct('rawAccel',true, 'scatterCorr',false, 'runningLikelihoodVsWinSize',true, 'runningLikelihoodFull',true));
		
		%% Play a video of the simulation
		if true, generateDronesInRoomVideo([], posUAVgt, roomDimensions, spotterCam); end
	end
end
