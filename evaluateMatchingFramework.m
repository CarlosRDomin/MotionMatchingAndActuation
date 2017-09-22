%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs a set of experiments with our matching framework and tries to determine which drone is which
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;
cdToThisScriptsDirectory;

simOutputFolder = [fileparts(mfilename('fullpath')) '/data/Simulations/'];
dims = 1:3;		% 3-axis acceleration
spotterCam = struct('fps',30, 'FOV',120*pi/180, 'pos',[], 'orient',[1 0 0; 0 0 -1; 0 1 0]);	% Camera's Field Of View in rad (120º), orientation pointing in the direction of the y-axis
threshRisk = 0.25; % m, how close drones can get to a wall (to avoid going through them or crashing)
threshPosteriorsEndsExperiment = 0.9999; % End the experiment when all drones are identified with confidence over 99.99%
maxSpeed = 1;	% m/s max speed drones can fly at
deltaT = 1;		% Timestep/iteration: 1s
deltaP = maxSpeed*deltaT;	% Max distance drones can move during a timestep
frameworkWinSize = deltaT*spotterCam.fps; iW=1;

sigmaNoiseCam = 0.10;	% m, noise in camera's position estimation
sigmaNoisePos = 0.05;	% m, error in position after performing a motion command
sigmaNoiseAccel = 0.25;	% m/s2, noise in IMU's accelerometer data
derivFiltOrder = 2; derivFiltHalfWinSize = 10;
typeOfExperiment = 'MatchingFramework';
repsPerExperiment = 20;

PLOT_EXPERIMENT = false;
SAVE_EXPERIMENT = true;

for roomDimensionsCell = {[15, 10, 3]} %{[5, 5, 2.5], [15, 10, 3]} % Width x Depth x Height of the room in m
	roomDimensions = roomDimensionsCell{:};
	spotterCam.pos = [roomDimensions(1)/2, -(roomDimensions(1)/2) / tan(spotterCam.FOV/2), roomDimensions(3)/2]; % Need FOV to determine pos

	for N = 2:50
		M = N;
		for normalizationByRowAndColumn = 0:2
			for typeOfMotionCell = {'random', 'hovering'} % {'random', 'hovering', 'landed'}
				typeOfMotion = typeOfMotionCell{:};

				% Initialize struct where we keep all fixed (constant) parameters
				paramFields = {'typeOfExperiment','typeOfMotion','N','normalizationByRowAndColumn','spotterCam','maxSpeed','deltaT','frameworkWinSize','dims','roomDimensions','threshRisk','sigmaNoiseCam','sigmaNoisePos','sigmaNoiseAccel','derivFiltOrder','derivFiltHalfWinSize'};
				paramStruct = cell2struct(cell(1, length(paramFields)), paramFields, 2);
				for f = paramFields	% Populate param struct
					paramStruct.(f{:}) = eval(f{:});
				end
				% Initialize struct where we'll save the experiment results
				variableFields = {'groundTruthAssignment','accelUAV','accelCam','posUAVgt','posUAVcam','t','runningWinScore','runningLikelihood','runningPosterior','assignedMatch','experimentRanOutOfTime'};
				variableStruct = cell2struct(cell(0, length(variableFields)), variableFields, 2);

				for experimentRep = 1:repsPerExperiment
					tMax = 300; t = 0 : 1/spotterCam.fps : tMax;
					experimentRanOutOfTime = true; % Initialize this variable to determine if an experiment concluded because all posteriors were above threshPosteriorsEndsExperiment or we ran out of time before that
					groundTruthAssignment = randperm(N,M)'; groundTruthAssignment=(1:M)'; % Get a random permutation of M elements picked from the set 1:N

					% Initialize random positions, making sure no one is too close to walls
% 					distToWalls = -1;
% 					while min(distToWalls(:)) < threshRisk % Keep trying if necessary until no one goes through walls
						if strcmp(typeOfMotion, 'landed')
							initPosUAV = cat(3, threshRisk + rand(N,1,length(dims)-1).*reshape(roomDimensions(1:length(dims)-1)-2*threshRisk, 1,1,[]), zeros(N,1,1));
						elseif strcmp(typeOfMotion, 'hovering')
							initPosUAV = cat(3, threshRisk + rand(N,1,length(dims)-1).*reshape(roomDimensions(1:length(dims)-1)-2*threshRisk, 1,1,[]), 0.5*roomDimensions(end)*ones(N,1,1));
						elseif strcmp(typeOfMotion, 'random')
							initPosUAV = threshRisk + rand(N,1,length(dims)).*reshape(roomDimensions-2*threshRisk, 1,1,[]);
						end
% 						[~,~,distToWalls] = estimateRiskOfCommand(zeros(3*N,1), [(1:N)' groundTruthAssignment], initPosUAV, roomDimensions);
% 					end
					posUAVgt = cat(2, initPosUAV, NaN(N, length(t)-1, length(dims)));
					posUAVcam = posUAVgt;

					% Init all other variables
					accelUAV = cat(2, zeros(N, 1, length(dims)), NaN(N, length(t)-1, length(dims)));
					accelCam = cat(2, zeros(M, derivFiltHalfWinSize, length(dims)), NaN(M, length(t)-derivFiltHalfWinSize, length(dims)));
					runningWinScore = cat(3, zeros(N, M, 1, length(dims), length(frameworkWinSize)), NaN(N, M, length(t)-1, length(dims), length(frameworkWinSize)));
					runningLikelihood = NaN(N, M, length(t), length(frameworkWinSize));
					runningPosterior = cat(3, ones(N, M, 1, length(frameworkWinSize))./(N+M-1), NaN(N, M, length(t)-1, length(frameworkWinSize)));
					assignedMatch = NaN(N, length(t), length(frameworkWinSize));

					if PLOT_EXPERIMENT
						figure('Units','normalized', 'Position',[0.3 0.4 0.4 0.25]);
						plotDronesInRoom(posUAVgt(:,1,:), roomDimensions, spotterCam);
					end

					for currT = 0 : deltaT : tMax-(1/spotterCam.fps)
						currTind = currT*spotterCam.fps + 1;

						% Check for a condition to end the experiment (everyone identified)
						if testIfExperimentShouldEnd(runningPosterior(:,:,currTind,iW), threshPosteriorsEndsExperiment)
							dispImproved(sprintf('\nAll drones identified >= %.2f%%! Concluding experiment at t=%.2fs\n', 100*threshPosteriorsEndsExperiment, currT), 'keepthis');
							accelUAV(:,currTind+1:end,:) = [];
							accelCam(:,currTind+1:end,:) = [];
							posUAVgt(:,currTind+1:end,:) = [];
							posUAVcam(:,currTind+1:end,:) = [];
							t(currTind+1:end) = [];
							runningWinScore(:,:,currTind+1:end,:,:) = [];
							runningLikelihood(:,:,currTind+1:end,:) = [];
							runningPosterior(:,:,currTind+1:end,:) = [];
							assignedMatch(:,currTind+1:end,:) = [];
							experimentRanOutOfTime = false;
							break;
						end

						posUAVcam(:,currTind,:) = posUAVgt(:,currTind,:) + sigmaNoiseCam.*randn(size(posUAVcam(:,currTind,:)));
						if strcmp(typeOfMotion, 'landed')
							a = zeros(N, frameworkWinSize, length(dims));
							p = zeros(N, frameworkWinSize, length(dims));
							posNoise = 0;
						elseif strcmp(typeOfMotion, 'hovering')
							posNoise = sigmaNoisePos;
							deltaToInitHoverPoint = squeeze(initPosUAV-posUAVcam(:,currTind,:));
							[az,el,rho] = cart2sph(deltaToInitHoverPoint(:,1), deltaToInitHoverPoint(:,2), deltaToInitHoverPoint(:,3));
							command = reshape([rho'; az'; pi/2-el'], [],1); % zeros(3*N,1);
							[a,~,p] = predictPosVelAccelFromCommand(command(1:3:end), command(2:3:end), command(3:3:end), spotterCam.fps, 1, deltaT); a=permute(a,3:-1:1); p=permute(p,3:-1:1);
						elseif strcmp(typeOfMotion, 'random')
							posNoise = sigmaNoisePos;
							distToWalls = -1;
							while min(distToWalls(:)) < threshRisk % Keep trying if necessary until no one goes through walls
								command = reshape([deltaP*rand(1,N); [2*pi; pi].*rand(2,N)], [],1); % Generate a random command
								[~,~,distToWalls] = estimateRiskOfCommand(command, [(1:N)' groundTruthAssignment], posUAVcam(:,currTind,:), roomDimensions);
							end
							[a,~,p] = predictPosVelAccelFromCommand(command(1:3:end), command(2:3:end), command(3:3:end), spotterCam.fps, 1, deltaT); a=permute(a,3:-1:1); p=permute(p,3:-1:1);
						end
						posUAVgt(:,currTind+(1:frameworkWinSize),:) = posUAVgt(:,currTind,:) + p + posNoise.*randn(N,frameworkWinSize,length(dims));
						posUAVcam(:,currTind+(1:frameworkWinSize),:) = posUAVgt(:,currTind+(1:frameworkWinSize),:) + sigmaNoiseCam.*randn(M,frameworkWinSize,length(dims));
						accelUAV(:,currTind+(1:frameworkWinSize),:) = a + sigmaNoiseAccel.*randn(N,frameworkWinSize,length(dims));
						%accelCam(:,currTind+(1:frameworkWinSize),:) = accelCam(:,currTind+(1:frameworkWinSize),:) + a(groundTruthAssignment,:,:) + sigmaNoiseAccel.*randn(M,frameworkWinSize,length(dims));
						tDerivFilterInds = currTind+(-2*derivFiltHalfWinSize:frameworkWinSize); tDerivFilterInds(tDerivFilterInds<1) = []; % Make sure this doesn't break on the first time window
						validDerivFilterInds = derivFiltHalfWinSize+1 : length(tDerivFilterInds)-derivFiltHalfWinSize; % First and last derivFiltHalfWinSize points are incorrect -> Don't overwrite them. Also, min(currTind, derivFiltHalfWinSize+1) forces to write indices 1:derivFiltHalfWinSize on the first time-window (currTind=1)
						filteredAccelCam = derivFilter(posUAVcam(:,tDerivFilterInds,:), 2, spotterCam.fps, derivFiltOrder, 2*derivFiltHalfWinSize+1);  % Save the results of the filter so we can index it from the derivFiltHalfWinSize+1'th point on
						accelCam(:,tDerivFilterInds(validDerivFilterInds),:) = filteredAccelCam(:,validDerivFilterInds,:);
						[runningWinScore, runningLikelihood, runningPosterior, assignedMatch] = computeBayesianIteration(runningWinScore, runningLikelihood, runningPosterior, assignedMatch, accelCam, accelUAV, currTind+frameworkWinSize, dims, frameworkWinSize, N, M, derivFiltHalfWinSize, [], normalizationByRowAndColumn);
						dispImproved(sprintf('\nPosterior likelihood: (@%s, motion="%s", N=%2d, norm=%d, experiment=%2d, currT=%3d)\n%s%s\n', datestr(now, 'HH:MM:SS'), typeOfMotion, N, normalizationByRowAndColumn, experimentRep, currT, [num2str(100.*runningPosterior(:,:,currTind+frameworkWinSize,iW), '%8.2f')'; repmat(13,1,N)], repmat('-',1,50)), 'keepthis');
						if PLOT_EXPERIMENT, plotDronesInRoom(posUAVgt(:,currTind+frameworkWinSize,:), roomDimensions, spotterCam); end
					end

					%% Save experiment repetition results in temp struct
					variableStruct(end+1).t = t; % Create a new row in the struct by setting any field
					for f = variableFields
						variableStruct(end).(f{:}) = eval(f{:});
					end
					% Just in case, save the *.mat (same name, so keep overwriting as we perform more repetitions)
					if SAVE_EXPERIMENT
						save_fileName = strjoin([simOutputFolder typeOfExperiment '_' typeOfMotion '_' N 'N_norm' normalizationByRowAndColumn '_' string(roomDimensions).join('x') '.mat'], '');
						save(save_fileName, 'paramStruct','variableStruct', '-v7.3'); % The flag '-v7.3' allows to save files of size >= 2GB
						dispImproved(sprintf('Saved simulation results as "%s"\n', save_fileName), 'keepthis');
					end
				end

				dispImproved(sprintf('\n\nFinished %d experiments with N=%d, typeOfMotion="%s"\n\n', repsPerExperiment, N, typeOfMotion), 'keepthis');

		% 		%% Plot last experiment results (raw accels, likelihood, posteriors...)
		% 		outFields = {'runningWinScore','runningLikelihood','runningPosterior','assignedMatch','N','M','iCams','dims','runningCorrWinSizes','t','accelCam','accelUAV'};
		% 		runningCorrStruct = cell2struct(cell(1, length(outFields)), outFields, 2);
		% 		iCams=groundTruthAssignment; runningCorrWinSizes=frameworkWinSize;
		% 		for f = outFields	% Populate output struct with results
		% 			runningCorrStruct.(f{:}) = eval(f{:});
		% 		end
		% 		plotExperimentResults(runningCorrStruct, struct('rawAccel',true, 'scatterCorr',false, 'runningLikelihoodVsWinSize',false, 'runningLikelihoodFull',true));
		% 		
		% 		%% Play a video of the simulation
		% 		generateDronesInRoomVideo([], posUAVgt, roomDimensions, spotterCam, false, frameworkWinSize);
			end
		end
	end
end
