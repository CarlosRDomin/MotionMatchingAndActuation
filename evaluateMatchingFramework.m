%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs a set of experiments with our matching framework and tries to determine which drone is which
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;
cdToThisScriptsDirectory;

dataOutputFolder = [fileparts(mfilename('fullpath')) '/data/'];
dims = 1:3;		% 3-axis acceleration
spotterCam = struct('fps',30, 'FOV',120*pi/180, 'pos',[], 'orient',[0 1 0; 0 0 -1; -1 0 0]);	% Camera's Field Of View in rad (120º), orientation pointing in the direction of the y-axis
threshRisk = 0.25; % m, how close drones can get to a wall (to avoid going through them or crashing)
threshPosteriorsEndsExperiment = 0.999; % End the experiment when all drones are identified with confidence over 99.99%
maxSpeed = 1;	% m/s max speed drones can fly at
deltaT = 1;		% Timestep/iteration: 1s
deltaP = maxSpeed*deltaT;	% Max distance drones can move during a timestep
frameworkWinSize = deltaT*spotterCam.fps; iW=1;

normalizationByRowAndColumn = 0; % Regular Bayesian normalization (divide each cell by the sum of its row)
sigmaNoiseCam = 0.05;			% m, noise in camera's position estimation
sigmaNoiseMotion = 0.05*deltaT;	% m/s2, error in each UAV's acceleration while performing a motion command
sigmaNoiseAccel = 0.25;			% m/s2, noise in IMU's accelerometer data
sigmaLikelihood = 1;
derivFiltOrder = 2; derivFiltHalfWinSize = 15;
typeOfExperiment = 'Simulation';
repsPerExperiment = 20;

PLOT_EXPERIMENT = false;
SAVE_EXPERIMENT = true;

for roomDimensionsScale = [1 2]
	roomDimensions = roomDimensionsScale * [2, 2, 1]; % roomDimensionsCell{:};
	spotterCam.pos = [roomDimensions(1) + (roomDimensions(1)/2) / tan(spotterCam.FOV/2), roomDimensions(2)/2, roomDimensions(3)/2]; % Need FOV to determine pos
	% Old way of positioning the camera (looking at y-axis): spotterCam.pos = [roomDimensions(1)/2, -(roomDimensions(1)/2) / tan(spotterCam.FOV/2), roomDimensions(3)/2]; % Need FOV to determine pos

	if roomDimensionsScale == 2
		Narray = ceil(24./(5:-1:1));
		Narray = ceil(24./(1:5));
	else
		Narray = 3;
	end
	for N = Narray
		M = N;
		for sigmaLikelihood = 1
			for typeOfMotionCell = {'oneAtATime', } %{'hovering', 'random', 'landed', 'oneAtATime', 'lowestRisky', 'oursSim'}
				typeOfMotion = typeOfMotionCell{:};
				if startsWith(typeOfMotion,'ours', 'IgnoreCase',true)
					sigmaNoiseMotion = 0.1*deltaT;	% m/s2, error in each UAV's acceleration while performing a motion command
					ourActuationNumLowRiskIterations = 5;		% Move all in the same dir, diff. amplitude during 5 iterations, then optimally
					ourActuationNumBestAssignments = max(3, ceil(N/2)); %max(1, min(4, ceil(N/3)));	% Number of most likely assignments to consider when generating command
					ourActuationStopLowRiskCriteria = 0.5;		% Stop lowRisk mode when P(assignment #ourActuationNumBestAssignments) < ourActuationStopLowRiskCriteria*P(assignment #1)
					paramsOurActuation = {'ourActuationNumLowRiskIterations', 'ourActuationNumBestAssignments', 'ourActuationStopLowRiskCriteria'};
					if strcmpi(typeOfMotion, 'oursReal')
						repsPerExperiment = 10;
						tMin = 20;
						tMax = 20;
						movingAvgFiltWinSize = 30;
					else
						tMin = 30;
						tMax = 30;
						movingAvgFiltWinSize = 2;
					end
				else
					sigmaNoiseMotion = 0.04*deltaT;	% m/s2, error in each UAV's acceleration while performing a motion command
					sigmaNoiseCam = 0.05;
					if strcmpi(typeOfMotion, 'hovering'), sigmaNoiseMotion = 0.01*deltaT; sigmaNoiseCam = 0.01; elseif strcmpi(typeOfMotion, 'random'), sigmaNoiseCam = 0.1; end
					threshPosteriorsEndsExperiment = 0.999; % End the experiment when all drones are identified with confidence over 99.99%
					tMin = 30;
					tMax = 150;
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
					if strcmpi(typeOfMotion, 'hovering') || false %|| strcmpi(typeOfMotion, 'lowestRisky')
						zScale = 0.5*roomDimensions(end);
					else % Any other combination, start landed
						zScale = 0;
					end
					if N>20, wallMargin = 0.8; else, wallMargin = 1; end % With high N's it's almost impossible to generate a valid combination, so use more space closer to the walls
					if strcmpi(typeOfMotion, 'random') || false
						initPosUAV = reshape(generateSpreadRandomPoints(N, threshRisk, roomDimensions, wallMargin), N,1,[]);
					else
						initPosUAV = cat(3, reshape(generateSpreadRandomPoints(N, threshRisk, roomDimensions(1:2), wallMargin), N,1,[]), zScale*ones(N,1,1));
					end
					posUAVgt = cat(2, initPosUAV, NaN(N, length(t)-1, length(dims)));
					posUAVcam = posUAVgt;
					velUAVgt = zeros(N, length(t), length(dims));
					accelUAVgt = zeros(N, length(t), length(dims));
					
					% Write the initial positions to a txt file for the Python script if evaluating our actuation in real life
					if strcmpi(typeOfMotion, 'oursReal')
						outputFolder = [dataOutputFolder 'Real/Ours/experiment' num2str(experimentRep) '/'];
						if exist(outputFolder, 'file') == 7
							error(['Couldn''t create folder' outputFolder ', does it already exist? Don''t want to overwrite data :)']);
						else
							mkdir(outputFolder);
						end
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
						% Create the folder for this iteration (if collecting real data for our actuation) so writing newCommand.txt doesn't fail
						if strcmpi(typeOfMotion, 'oursReal')
							iterOutputFolder = [outputFolder 'iteration' num2str(currIter) '/'];
							if exist(iterOutputFolder, 'file') == 7
								error(['Couldnt create folder' iterOutputFolder]);
							else
								mkdir(iterOutputFolder)
							end
						end

						% Check for a condition to end the experiment (everyone identified)
						if currT > tMin && testIfExperimentShouldEnd(runningPosterior(:,:,currTind,iW), threshPosteriorsEndsExperiment)
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
							p = zeros(N, frameworkWinSize, length(dims)) + posUAVgt(:,currTind,:);
						else
							noiseMotion = sigmaNoiseMotion;
							if strcmpi(typeOfMotion, 'hovering') || (strcmpi(typeOfMotion, 'oneAtATime') && currIter > N)
								desiredPosUAV = initPosUAV; desiredPosUAV(:,:,end) = roomDimensions(end)/2;
								% if strcmpi(typeOfMotion, 'oneAtATime'), desiredPosUAV(:,:,end) = roomDimensions(end)/2; end
								deltaToInitHoverPoint = squeeze(initPosUAV-posUAVcam(:,currTind,:));
								[az,el,rho] = cart2sph(deltaToInitHoverPoint(:,1), deltaToInitHoverPoint(:,2), deltaToInitHoverPoint(:,3));
								command = reshape([rho'; az'; pi/2-el'], [],1); % zeros(3*N,1);
							elseif strcmpi(typeOfMotion, 'random')
								if currTind == 1
									maxRho = min(roomDimensions(3)/2, deltaP); dirTheta = 0; dirPhi = 0; % Move up
									command = reshape(repmat([maxRho; dirTheta; dirPhi], 1,N), [],1);
								else
									command = reshape([deltaP*rand(1,N); [2*pi; pi].*rand(2,N)], [],1); % Generate a random command
% 									distToWalls = -1;
% 									while min(distToWalls(:)) < threshRisk % Keep trying if necessary until no one goes through walls
% 										command = reshape([deltaP*rand(1,N); [2*pi; pi].*rand(2,N)], [],1); % Generate a random command
% 										[~,~,distToWalls] = estimateRiskOfCommand(command, [(1:N)' groundTruthAssignment], posUAVcam(:,currTind,:), roomDimensions);
% 									end
								end
							elseif strcmpi(typeOfMotion, 'lowestRisky')
								[rhoBounds, dirTheta, dirPhi] = estimateLowestRiskyDirOfMotion(reshape(posUAVcam(:,currTind,:), [],size(posUAVcam,3)), roomDimensions, threshRisk, deltaP);
								command = reshape([sum(rhoBounds)*rand(1,N)-rhoBounds(2); repmat([dirTheta; dirPhi], 1,N)], [],1);	% All drones move with same dirTheta and dirPhi, random rhos (of up to maxRho)
								%command = reshape([maxRho*((randperm(N)-1)/(N-1)); repmat([dirTheta; dirPhi], 1,N)], [],1);	% All drones move with same dirTheta and dirPhi, random rhos (of up to maxRho)
							elseif strcmpi(typeOfMotion, 'oneAtATime')
								maxRho = min(roomDimensions(3)/2, deltaP); dirTheta = 0; dirPhi = 0; % Move up
								command = reshape([maxRho*((1:N)==mod(currT/deltaT, N)+1); repmat([dirTheta; dirPhi], 1,N)], [],1);
							elseif startsWith(typeOfMotion,'ours', 'IgnoreCase',true)
								runningPriorTemp = runningPosterior(:,:,currTind,iW);
								[rhoBounds, dirTheta, dirPhi] = estimateLowestRiskyDirOfMotion(reshape(posUAVcam(:,currTind,:), [],size(posUAVcam,3)), roomDimensions, threshRisk, deltaP);
								
								%sortedPriors = sort(runningPriorTemp,2);
								%if all(sortedPriors(:,end-1) > ourActuationStopLowRiskCriteria*sortedPriors(:,end))
								
								assignmentList = computeNBestAssignments(ourActuationNumBestAssignments, runningPriorTemp, -log(1e-30),-log(1e-2));
								%P_assignment_1 = prod(runningPriorTemp(sub2ind(size(runningPriorTemp), assignmentList(1).matches(:,1),assignmentList(1).matches(:,2))));
								%P_assignment_k = prod(runningPriorTemp(sub2ind(size(runningPriorTemp), assignmentList(ourActuationNumBestAssignments).matches(:,1),assignmentList(ourActuationNumBestAssignments).matches(:,2))));
								%if P_assignment_k > ourActuationStopLowRiskCriteria*P_assignment_1 % Keep doing low risk until k-th assignment has lower prob than 0.1*top assignment
								if currIter <= ourActuationNumLowRiskIterations
									command = reshape([sum(rhoBounds)*rand(1,N)-rhoBounds(2); repmat([dirTheta; dirPhi], 1,N)], [],1);	% All drones move with same dirTheta and dirPhi, random rhos (of up to maxRho)
								else
									P_total = 0;
									expectedCommand = zeros(N,3);
									startingCommand = commandWithValidBounds(reshape([sum(rhoBounds)*rand(1,N)-rhoBounds(2); repmat([dirTheta; dirPhi], 1,N)], [],1));
									
									for i = 1:ourActuationNumBestAssignments
										assignments = assignmentList(i).matches; unassignedUAVs = assignmentList(i).unassignedUAVs; unassignedDetections = assignmentList(i).unassignedDetections;
										P_assignment = prod(runningPriorTemp(sub2ind(size(runningPriorTemp), assignments(:,1),assignments(:,2))));
										if P_assignment<0.2*P_total/(ourActuationNumBestAssignments-1), break; end % Stop if this combination is already unlikely
										P_total = P_total + P_assignment;
										%sortedAssignment = sortrows([assignments; [unassignedUAVs(:) NaN(length(unassignedUAVs),1)]]);
										[command,newP,exitFlag,output] = findOptimalCommandGivenAssignment(assignments, posUAVcam(:,currTind,:), roomDimensions, threshRisk, runningPriorTemp, spotterCam, deltaP, deltaT, sigmaNoiseAccel, startingCommand);
										expectedCommand = expectedCommand + P_assignment.*reshape(command, 3,[])';
									end
									command = deltaPtoCommand(expectedCommand./P_total);
								end
							end
							[a,v,p] = predictPosVelAccelFromCommand(command(1:3:end), command(2:3:end), command(3:3:end), accelUAVgt(:,currTind,:), velUAVgt(:,currTind,:), posUAVgt(:,currTind,:), spotterCam.fps, deltaT, noiseMotion);
						end
						
						% Update pos and accel based on the new command (or wait to receive data from Python if it's a real experiment)
						if strcmpi(typeOfMotion, 'oursReal')
							% First, write the command to a file so Python can tell each drone where to move
							fileID = fopen([iterOutputFolder 'newCommand.txt'], 'w');
							commandDeltaP = commandToDeltaP(command(1:3:end), command(2:3:end), command(3:3:end)); % N x 1 x 3
							strDeltaP = num2str(commandDeltaP, '%.3f,');  % Convert commandDeltaP to string and save it in a temporary variable so we can remove the trailing ','
							strCurrentPos = num2str(posUAVgt(:,currTind,:), '%.3f,');  % Convert current positions to string
							fprintf(fileID, strcat(strCurrentPos, strDeltaP(:,1:end-1), '\n')');  % Separate rows with a '\n' and transpose (otherwise it would write column-wise)
							fprintf(fileID, 'Done');  % Make sure we signal Python we're done writing the file (avoids race conditions where Python only reads when we're half way through writing)
							fclose(fileID);
							
							% Wait for Python to collect data and write log files
							dispImproved('New command issued, waiting for Python...');
							while ~exist([iterOutputFolder 'dataCollectionDone.txt'], 'file')
								pause(1);	% Keep polling every second
							end
							
							% Load logs
							while true
								try
									data = loadRealExperimentData(struct('datetime',{['iteration' num2str(currIter)]}, 'ch',['experiment' num2str(experimentRep)]), outputFolder(1:end-1), derivFiltOrder, 1+2*derivFiltHalfWinSize, movingAvgFiltWinSize);
									if isempty(find(data.drone_id.measured == N, 1)) % Check that the data was correctly recorded (sometimes there's errors before the iteration ends)
										dispImproved('ERROR, wrong data collected, recollect this iteration');
										pause(5);
									else % All good :)
										break;
									end
								catch
									pause(5);
								end
							end
							
							for n = 1:N
								ind_start = find(data.drone_id.measured == n, 1);
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
						save_fileName = strjoin([outputFolder typeOfExperiment '_' typeOfMotion '_' N 'N_' string(roomDimensions).join('x') '.mat'], '');
						save(save_fileName, 'paramStruct','variableStruct', '-v7.3'); % The flag '-v7.3' allows to save files of size >= 2GB
						dispImproved(sprintf('Saved simulation results as "%s"\n', save_fileName), 'keepthis');
					end

					if false
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
						
						%% Compute survival rate
						threshCollision=0.25;
						distAmongDrones = zeros(N*(N-1)/2, currTind);
						survivingDrones = true(N, currTind);
						avoidDiagDist = threshCollision.*eye(N); % Use this helper matrix to avoid finding drone crashes of a drone with itself
						prevSurvivors = true(N,1);
						for tInd = 1:currTind
							distAmongDrones(:,tInd) = pdist(reshape(posUAVgt(:,tInd,:), [],size(posUAVgt,3)));
							[deadDronesI, deadDronesJ] = find((squareform(distAmongDrones(:,tInd)) + avoidDiagDist) < threshCollision);
							shouldDie = unique([deadDronesI; deadDronesJ]);
							shouldDie(prevSurvivors(shouldDie)==false) = []; % Dead drones can't "kill" other drones
							prevSurvivors(shouldDie) = 0;
							survivingDrones(:,tInd) = prevSurvivors;
						end
						survivalRate = sum(survivingDrones,1)/N;
						
						figure;
						subplot(2,1,1); plot(min(distAmongDrones,[],1)); title('Minimum drone-drone distance');
						subplot(2,1,2); plot(survivalRate); title('Survival rate');
					end
				end
				dispImproved(sprintf('\n\nFinished %d experiments with N=%d, typeOfMotion="%s"\n\n', repsPerExperiment, N, typeOfMotion), 'keepthis');
			end
		end
	end
end
