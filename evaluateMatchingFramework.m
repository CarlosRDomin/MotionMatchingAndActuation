%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs a set of experiments with our matching framework and tries to determine which drone is which
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;
cdToThisScriptsDirectory;

simOutputFolder = [fileparts(mfilename('fullpath')) '/data/Simulations/'];
dims = 1:3;		% 3-axis acceleration
spotterCam = struct('fps',30, 'FOV',120*pi/180, 'pos',[], 'orient',[1 0 0; 0 0 -1; 0 1 0]);	% Camera's Field Of View in rad (120º), orientation pointing in the direction of the y-axis
threshRisk = 0.5; % m, how close drones can get to a wall (to avoid going through them or crashing)
threshPosteriorsEndsExperiment = 0.9999; % End the experiment when all drones are identified with confidence over 99.99%
maxSpeed = 1;	% m/s max speed drones can fly at
deltaT = 1;		% Timestep/iteration: 1s
deltaP = maxSpeed*deltaT;	% Max distance drones can move during a timestep
frameworkWinSize = deltaT*spotterCam.fps; iW=1;

normalizationByRowAndColumn = 0; % Regular Bayesian normalization (divide each cell by the sum of its row)
sigmaNoiseCam = 0.05;	% m, noise in camera's position estimation
sigmaNoisePos = 0.01;	% m, error in position after performing a motion command
sigmaNoiseAccel = 0.25;	% m/s2, noise in IMU's accelerometer data
sigmaLikelihood = 1;
derivFiltOrder = 2; derivFiltHalfWinSize = 12;
typeOfExperiment = 'MatchingFramework';
repsPerExperiment = 20;

PLOT_EXPERIMENT = false;
SAVE_EXPERIMENT = false;

for roomDimensionsCell = {[5, 5, 2.5]} %{[5, 5, 2.5], [15, 10, 3]} % Width x Depth x Height of the room in m
	roomDimensions = roomDimensionsCell{:};
	spotterCam.pos = [roomDimensions(1)/2, -(roomDimensions(1)/2) / tan(spotterCam.FOV/2), roomDimensions(3)/2]; % Need FOV to determine pos

	for N = 5
		M = N;
		for sigmaLikelihood = 1 %[0.05:0.05:0.2 0.25:0.25:1.5]
			for typeOfMotionCell = {'oneAtATime', 'hovering'} %{'random', 'hovering', 'landed', 'oneAtATime', 'lowestRisky', 'ours'}
				typeOfMotion = typeOfMotionCell{:};

				% Initialize struct where we keep all fixed (constant) parameters
				paramFields = {'typeOfExperiment','typeOfMotion','N','normalizationByRowAndColumn','spotterCam','maxSpeed','deltaT','frameworkWinSize','dims','roomDimensions','threshRisk','sigmaNoiseCam','sigmaNoisePos','sigmaNoiseAccel','sigmaLikelihood','derivFiltOrder','derivFiltHalfWinSize'};
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
 					distAmongDrones = -1;
 					while min(distAmongDrones(:)) < threshRisk % Keep trying if necessary until no one goes through walls
						if strcmpi(typeOfMotion, 'hovering') %|| strcmpi(typeOfMotion, 'lowestRisky')
							initPosUAV = cat(3, threshRisk + rand(N,1,length(dims)-1).*reshape(roomDimensions(1:length(dims)-1)-2*threshRisk, 1,1,[]), 0.5*roomDimensions(end)*ones(N,1,1));
						elseif strcmpi(typeOfMotion, 'random')
							initPosUAV = threshRisk + rand(N,1,length(dims)).*reshape(roomDimensions-2*threshRisk, 1,1,[]);
						else % Any other combination, start landed
							initPosUAV = cat(3, threshRisk + rand(N,1,length(dims)-1).*reshape(roomDimensions(1:length(dims)-1)-2*threshRisk, 1,1,[]), zeros(N,1,1));
						end
 						[~,distAmongDrones,~] = estimateRiskOfCommand(zeros(3*N,1), [(1:N)' groundTruthAssignment], initPosUAV, roomDimensions);
 					end
					posUAVgt = cat(2, initPosUAV, NaN(N, length(t)-1, length(dims)));
					posUAVcam = posUAVgt;

					% Init all other variables
					accelUAV = cat(2, zeros(N, 1, length(dims)), NaN(N, length(t)-1, length(dims)));
					accelCam = cat(2, zeros(M, derivFiltHalfWinSize, length(dims)), NaN(M, length(t)-derivFiltHalfWinSize, length(dims)));
					runningWinScore = cat(3, zeros(N, M, 1, length(dims), length(frameworkWinSize)), NaN(N, M, length(t)-1, length(dims), length(frameworkWinSize)));
					runningLikelihood = NaN(N, M, length(t), length(frameworkWinSize));
					if normalizationByRowAndColumn==1, normalizeBy = N+M+1; else, normalizeBy = N; end
					runningPosterior = cat(3, ones(N, M, 1, length(frameworkWinSize))./normalizeBy, NaN(N, M, length(t)-1, length(frameworkWinSize)));
					assignedMatch = cat(2, repmat(randperm(N)', 1,1,length(frameworkWinSize)), NaN(N, length(t)-1, length(frameworkWinSize)));

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

						% posUAVcam(:,currTind,:) = posUAVgt(:,currTind,:) + sigmaNoiseCam.*randn(size(posUAVcam(:,currTind,:)));
						if strcmpi(typeOfMotion, 'landed')
							posNoise = 0;
							a = zeros(N, frameworkWinSize, length(dims));
							p = zeros(N, frameworkWinSize, length(dims));
						else
							posNoise = sigmaNoisePos;
							if strcmpi(typeOfMotion, 'hovering')
								deltaToInitHoverPoint = squeeze(initPosUAV-posUAVcam(:,currTind,:));
								[az,el,rho] = cart2sph(deltaToInitHoverPoint(:,1), deltaToInitHoverPoint(:,2), deltaToInitHoverPoint(:,3));
								command = reshape([rho'; az'; pi/2-el'], [],1); % zeros(3*N,1);
							elseif strcmpi(typeOfMotion, 'random')
								distToWalls = -1;
								while min(distToWalls(:)) < threshRisk % Keep trying if necessary until no one goes through walls
									command = reshape([deltaP*rand(1,N); [2*pi; pi].*rand(2,N)], [],1); % Generate a random command
									[~,~,distToWalls] = estimateRiskOfCommand(command, [(1:N)' groundTruthAssignment], posUAVcam(:,currTind,:), roomDimensions);
								end
							elseif strcmpi(typeOfMotion, 'lowestRisky')
								[maxRho, dirTheta, dirPhi] = estimateLowestRiskyDirOfMotion(reshape(posUAVcam(:,currTind,:), [],size(posUAVcam,3)), roomDimensions, threshRisk);
								maxRho = min(maxRho, deltaP); % Make sure it's a feasible command
								command = reshape([maxRho*rand(1,N); repmat([dirTheta; dirPhi], 1,N)], [],1);	% All drones move with same dirTheta and dirPhi, random rhos (of up to maxRho)
								%command = reshape([maxRho*((randperm(N)-1)/(N-1)); repmat([dirTheta; dirPhi], 1,N)], [],1);	% All drones move with same dirTheta and dirPhi, random rhos (of up to maxRho)
							elseif strcmpi(typeOfMotion, 'oneAtATime')
								[maxRho, dirTheta, dirPhi] = estimateLowestRiskyDirOfMotion(reshape(posUAVcam(:,currTind,:), [],size(posUAVcam,3)), roomDimensions, threshRisk);
								maxRho = min(maxRho, deltaP); % Make sure it's a feasible command
								command = reshape([maxRho*((1:N)==mod(currT/deltaT, N)+1); repmat([dirTheta; dirPhi], 1,N)], [],1);
							end
							[a,~,p] = predictPosVelAccelFromCommand(command(1:3:end), command(2:3:end), command(3:3:end), spotterCam.fps, 1, deltaT); a=permute(a,3:-1:1); p=permute(p,3:-1:1);
						end
						posUAVgt(:,currTind+(1:frameworkWinSize),:) = posUAVgt(:,currTind,:) + p + posNoise.*randn(N,frameworkWinSize,length(dims));
						posUAVcam(:,currTind+(1:frameworkWinSize),:) = posUAVgt(:,currTind+(1:frameworkWinSize),:) + sigmaNoiseCam.*randn(M,frameworkWinSize,length(dims));
						accelUAV(:,currTind+(1:frameworkWinSize),:) = a + sigmaNoiseAccel.*randn(N,frameworkWinSize,length(dims));
						accelCam = estimateAccelCamFromPosCam(posUAVcam, accelCam, currTind, derivFiltOrder, derivFiltHalfWinSize, frameworkWinSize, spotterCam.fps);
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

					if true
						%% Plot last experiment results (raw accels, likelihood, posteriors...)
						outFields = {'runningWinScore','runningLikelihood','runningPosterior','assignedMatch','N','M','iCams','dims','runningCorrWinSizes','t','accelCam','accelUAV'};
						runningCorrStruct = cell2struct(cell(1, length(outFields)), outFields, 2);
						iCams=groundTruthAssignment; runningCorrWinSizes=frameworkWinSize;
						for f = outFields	% Populate output struct with results
							runningCorrStruct.(f{:}) = eval(f{:});
						end
						plotExperimentResults(runningCorrStruct, struct('rawAccel',true, 'scatterCorr',false, 'runningLikelihoodVsWinSize',false, 'runningLikelihoodFull',true));

						%% Play a video of the simulation
						generateDronesInRoomVideo([], posUAVgt, roomDimensions, spotterCam, false);
					end
				end
				dispImproved(sprintf('\n\nFinished %d experiments with N=%d, typeOfMotion="%s"\n\n', repsPerExperiment, N, typeOfMotion), 'keepthis');
			end
		end
	end
end
