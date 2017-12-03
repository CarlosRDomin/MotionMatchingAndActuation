%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs a set of experiments with our matching framework and tries to determine which drone is which
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;
cdToThisScriptsDirectory;

dataOutputFolder = [fileparts(mfilename('fullpath')) '/data/'];
dims = 1:3;		% 3-axis acceleration
spotterCam = struct('fps',30, 'FOV',120*pi/180, 'pos',[], 'orient',[1 0 0; 0 0 -1; 0 1 0]);	% Camera's Field Of View in rad (120º), orientation pointing in the direction of the y-axis
threshRisk = 0.5; % m, how close drones can get to a wall (to avoid going through them or crashing)
threshPosteriorsEndsExperiment = 0.99; % End the experiment when all drones are identified with confidence over 99.99%
maxSpeed = 1;	% m/s max speed drones can fly at
deltaT = 1;		% Timestep/iteration: 1s
deltaP = maxSpeed*deltaT;	% Max distance drones can move during a timestep
frameworkWinSize = deltaT*spotterCam.fps; iW=1;

normalizationByRowAndColumn = 0; % Regular Bayesian normalization (divide each cell by the sum of its row)
sigmaNoiseCam = 0.05;			% m, noise in camera's position estimation
sigmaNoiseMotion = 0.1*deltaT;	% m/s2, error in each UAV's acceleration while performing a motion command
sigmaNoiseAccel = 0.25;			% m/s2, noise in IMU's accelerometer data
sigmaLikelihood = 1;
derivFiltOrder = 2; derivFiltHalfWinSize = 15;
typeOfExperiment = 'Actuation';
repsPerExperiment = 20;

PLOT_EXPERIMENT = false;
SAVE_EXPERIMENT = false;

for roomDimensionsCell = {[8, 5, 3]} %{[5, 5, 2.5], [15, 10, 3]} % Width x Depth x Height of the room in m
	roomDimensions = roomDimensionsCell{:};
	spotterCam.pos = [roomDimensions(1)/2, -(roomDimensions(1)/2) / tan(spotterCam.FOV/2), roomDimensions(3)/2]; % Need FOV to determine pos

	for N = 10
		M = N;
		for sigmaLikelihood = 1 %[0.05:0.05:0.2 0.25:0.25:1.5]
			for typeOfMotionCell = {'oursReal', } %'oursSim', 'random', 'hovering', 'oneAtATime', 'lowestRisky', } %{'hovering', 'landed', 'oneAtATime', 'lowestRisky', 'ours'}
				typeOfMotion = typeOfMotionCell{:};
				if strcmpi(typeOfMotion, 'oursReal')
					tMax = 30;
					movingAvgFiltWinSize = 30;	% For IMU real data
					ourActuationNumLowRiskIterations = 5;		% Move all in the same dir, diff. amplitude during 5 iterations, then optimally
					ourActuationNumBestAssignments = ceil(N);	% Number of most likely assignments to consider when generating command
					paramsOurActuation = {'ourActuationNumLowRiskIterations', 'ourActuationNumBestAssignments'};
				else
					tMax = 300;
					movingAvgFiltWinSize = 2;	% For IMU simulated data
					paramsOurActuation = {};  % Don't save any parameters related to actuation
				end

				% Initialize struct where we keep all fixed (constant) parameters
				paramFields = [{'typeOfExperiment','typeOfMotion','N','normalizationByRowAndColumn','spotterCam','maxSpeed','deltaT','frameworkWinSize','dims','roomDimensions','threshRisk','sigmaNoiseCam','sigmaNoiseMotion','sigmaNoiseAccel','sigmaLikelihood','derivFiltOrder','derivFiltHalfWinSize','movingAvgFiltWinSize'}, paramsOurActuation];
				paramStruct = cell2struct(cell(1, length(paramFields)), paramFields, 2);
				for f = paramFields	% Populate param struct
					paramStruct.(f{:}) = eval(f{:});
				end
				% Initialize struct where we'll save the experiment results
				variableFields = {'groundTruthAssignment','accelUAV','accelCam','posUAVgt','posUAVcam','t','runningWinScore','runningLikelihood','runningPosterior','assignedMatch','experimentRanOutOfTime'};
				variableStruct = cell2struct(cell(0, length(variableFields)), variableFields, 2);

				for experimentRep = 1:repsPerExperiment
					t = 0 : 1/spotterCam.fps : tMax;
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
					velUAVgt = zeros(N, length(t), length(dims));
					accelUAVgt = zeros(N, length(t), length(dims));
					
					% Write the initial positions to a txt file for the Python script if evaluating our actuation in real life
					if strcmpi(typeOfMotion, 'oursReal')
						outputFolder = [dataOutputFolder 'Real/Ours/Experiment_' experimentRep '/'];
						fileID = fopen([outputFolder 'initPos.txt'], 'w');
						strInitPos = num2str(initPosUAV, '%.3f,');  % Convert initPosUAV to string and save it in a temporary variable so we can remove the trailing ','
						fprintf(fileID, strcat(strInitPos(:,1:end-1), '\n')');  % Separate rows with a '\n' and transpose (otherwise it would write column-wise)
						fclose(fileID);
					else
						outputFolder = [dataOutputFolder 'Simulations/'];
					end

					% Init all other variables
					[accelUAV, accelCam, runningWinScore, runningLikelihood, runningPosterior, assignedMatch] = initFrameworkVars(N, M, dims, t, frameworkWinSize, derivFiltHalfWinSize, normalizationByRowAndColumn);

					if PLOT_EXPERIMENT
						figure('Units','normalized', 'Position',[0.3 0.4 0.4 0.25]);
						plotDronesInRoom(initPosUAV, roomDimensions, spotterCam);
					end

					for currT = 0 : deltaT : tMax-(1/spotterCam.fps)
						currTind = currT*spotterCam.fps + 1;
						currIter = (currT/deltaT) + 1;
						iterOutputFolder = [outputFolder 'iteration_' num2str(currIter) '/'];

						% Check for a condition to end the experiment (everyone identified)
						if testIfExperimentShouldEnd(runningPosterior(:,:,currTind,iW), threshPosteriorsEndsExperiment)
							dispImproved(sprintf('\nAll drones identified >= %.2f%%! Concluding experiment at t=%.2fs\n', 100*threshPosteriorsEndsExperiment, currT), 'keepthis');
							accelUAV(:,currTind+1:end,:) = [];
							accelCam(:,currTind+1:end,:) = [];
							posUAVgt(:,currTind+1:end,:) = [];
							posUAVcam(:,currTind+1:end,:) = [];
							velUAVgt(:,currTind+1:end,:) = [];
							accelUAVgt(:,currTind+1:end,:) = [];
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
							noiseMotion = 0;
							a = zeros(N, frameworkWinSize, length(dims));
							v = zeros(N, frameworkWinSize, length(dims));
							p = zeros(N, frameworkWinSize, length(dims));
						else
							noiseMotion = sigmaNoiseMotion;
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
							elseif startsWith(typeOfMotion,'ours', 'IgnoreCase',true)
								runningPriorTemp = runningPosterior(:,:,currTind,iW);
								% if ~all(runningPriorTemp(sub2ind(size(runningPriorTemp), (1:N)',assignedMatch(:,currTind,iW))) > 1.5/N)
								if currIter <= ourActuationNumLowRiskIterations
									[maxRho, dirTheta, dirPhi] = estimateLowestRiskyDirOfMotion(reshape(posUAVcam(:,currTind,:), [],size(posUAVcam,3)), roomDimensions, threshRisk);
									maxRho = min(maxRho, deltaP); % Make sure it's a feasible command
									command = reshape([maxRho*rand(1,N); repmat([dirTheta; dirPhi], 1,N)], [],1);	% All drones move with same dirTheta and dirPhi, random rhos (of up to maxRho)
								else
									P_total = 0;
									expectedCommand = zeros(3*N,1);
									assignmentList = computeNBestAssignments(ourActuationNumBestAssignments, runningPriorTemp, -log(1e-30),-log(1e-2));

									for i = 1:ourActuationNumBestAssignments
										assignments = assignmentList(i).matches; unassignedUAVs = assignmentList(i).unassignedUAVs; unassignedDetections = assignmentList(i).unassignedDetections;
										P_assignment = prod(runningPriorTemp(sub2ind(size(runningPriorTemp), assignments(:,1),assignments(:,2))));
										if P_assignment<0.1*P_total/(ourActuationNumBestAssignments-1), break; end % Stop if this combination is already unlikely
										P_total = P_total + P_assignment;
										%sortedAssignment = sortrows([assignments; [unassignedUAVs(:) NaN(length(unassignedUAVs),1)]]);
										[command,newP,exitFlag,output] = fmincon(@(x) (-estimateImprovementOfCommand(x,assignments,runningPriorTemp,[],sigmaNoiseAccel,[],spotterCam.fps,deltaT)), ...
											zeros(3*N,1),[],[],[],[],zeros(3*N,1),repmat([deltaP,2*pi,pi]',N,1), ...
											@(x) deal(threshRisk-estimateRiskOfCommand(x,assignments,posUAVcam(:,currTind,:), roomDimensions), 0)); %, optimoptions('fmincon', 'Algorithm','active-set'));
										expectedCommand = expectedCommand + P_assignment.*command;
									end
									command = expectedCommand./P_total;
								end
							end
							[a,v,p] = predictPosVelAccelFromCommand(command(1:3:end), command(2:3:end), command(3:3:end), accelUAVgt(:,currTind,:), velUAVgt(:,currTind,:), posUAVgt(:,currTind,:), spotterCam.fps, deltaT, noiseMotion);
						end
						
						% Update pos and accel based on the new command (or wait to receive data from Python if it's a real experiment)
						if strcmpi(typeOfMotion, 'oursReal')
							% TODO
							% First, write the command to a file so Python can tell each drone where to move
							fileID = fopen([iterOutputFolder 'newCommand.txt'], 'w');
							commandDeltaP = commandToDeltaP(command(1:3:end), command(2:3:end), command(3:3:end)); % N x 1 x 3
							strDeltaP = num2str(commandDeltaP, '%.3f,');  % Convert commandDeltaP to string and save it in a temporary variable so we can remove the trailing ','
							strCurrentPos = num2str(posUAVgt(:,currTind,:), '%.3f,');  % Convert current positions to string
							fprintf(fileID, strcat(strCurrentPos, strDeltaP(:,1:end-1), '\n')');  % Separate rows with a '\n' and transpose (otherwise it would write column-wise)
							fclose(fileID);
							
							% Wait for Python to collect data and write log files
							while ~exist([iterOutputFolder 'dataCollectionDone.txt'], 'file')
								pause(1);	% Keep polling every second
							end
							
							% Load logs
							data = loadRealExperimentData(struct('datetime',{['iteration_' num2str(currIter)]}, 'ch',['exp' num2str(experimentRep)]), outputFolder, derivFiltOrder, 1+2*derivFiltHalfWinSize, movingAvgFiltWinSize);
							for n = 1:N
								ind_start = find(data.drone_id == n, 1);
								ind_drone_data = ind_start:(ind_start + frameworkWinSize - 1);
								for d = 1:size(accelUAV,3)
									accelUAV(n,currTind+(1:frameworkWinSize),d) = data.a_UAV.(char('X'+d-1)).measured(ind_drone_data);
									% Take camera data as ground truth (will be copied as posUAVcam and accelCam as well outside the for loop)
									accelUAVgt(n,currTind+(1:frameworkWinSize),d) = data.a_cam.(char('X'+d-1)).measured(ind_drone_data);
									velUAVgt(n,currTind+(1:frameworkWinSize),d) = data.v_cam.(char('X'+d-1)).measured(ind_drone_data);
									posUAVgt(n,currTind+(1:frameworkWinSize),d) = data.p_cam.(char('X'+d-1)).measured(ind_drone_data);
								end
							end
							posUAVcam(:,currTind+(1:frameworkWinSize),:) = posUAVgt(:,currTind+(1:frameworkWinSize),:);
							accelCam(:,currTind+(1:frameworkWinSize),:) = accelUAVgt(:,currTind+(1:frameworkWinSize),:);
						else
							accelUAVgt(:,currTind+(1:frameworkWinSize),:) = a;
							velUAVgt(:,currTind+(1:frameworkWinSize),:) = v;
							posUAVgt(:,currTind+(1:frameworkWinSize),:) = p;
							posUAVcam(:,currTind+(1:frameworkWinSize),:) = posUAVgt(:,currTind+(1:frameworkWinSize),:) + sigmaNoiseCam.*randn(M,frameworkWinSize,length(dims));
							accelUAV(:,currTind+(1:frameworkWinSize),:) = a + sigmaNoiseAccel.*randn(N,frameworkWinSize,length(dims));
							%for n=1:size(accelUAV,1), for d=1:size(accelUAV,3), accelUAV(n,currTind+(1:frameworkWinSize),d) = sgolayfilt(accelUAV(n,currTind+(1:frameworkWinSize),d), 1, 5); end; end
							for n=1:size(accelUAV,1), for d=1:size(accelUAV,3), accelUAV(n,currTind+(1:frameworkWinSize),d) = movingAvgFilter(movingAvgFiltWinSize, accelUAV(n,currTind+(1:frameworkWinSize),d)); end; end
							accelCam = estimateAccelCamFromPosCam(posUAVcam, accelCam, currTind, derivFiltOrder, derivFiltHalfWinSize, frameworkWinSize, spotterCam.fps);
						end
						
						% Perform matching iteration
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
					if SAVE_EXPERIMENT || strcmpi(typeOfMotion, 'oursReal')  % Make sure we always save the data on real actuation!! Don't want to have to repeat the experiments...
						save_fileName = strjoin([outputFolder typeOfExperiment '_' typeOfMotion '_' N 'N_norm' normalizationByRowAndColumn '_' string(roomDimensions).join('x') '.mat'], '');
						save(save_fileName, 'paramStruct','variableStruct', '-v7.3'); % The flag '-v7.3' allows to save files of size >= 2GB
						dispImproved(sprintf('Saved simulation results as "%s"\n', save_fileName), 'keepthis');
					end

					if true
						%% Plot last experiment results (raw accels, likelihood, posteriors...)
						outFields = {'runningWinScore','runningLikelihood','runningPosterior','assignedMatch','N','M','iCams','dims','frameworkWinSize','t','accelCam','accelUAV'};
						runningCorrStruct = cell2struct(cell(1, length(outFields)), outFields, 2);
						iCams=groundTruthAssignment;
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
