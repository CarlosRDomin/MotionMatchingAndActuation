%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots all figures related to the matching framework evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;
cdToThisScriptsDirectory;

realOutputFolder = [fileparts(mfilename('fullpath')) '/data/Real/'];
simOutputFolder = [fileparts(mfilename('fullpath')) '/data/Simulations/'];
realOutputFolderAndPrefix = [realOutputFolder 'Real'];
simOutputFolderAndPrefix = [simOutputFolder 'Simulation'];
skipPreProcess = true; % Whether or not to "generate" by "overlapping" real experiment data
skipProcess = false; % Whether or not to (re)process all the experiment data

% Constants
confidenceThresholds = 0:0.025:1;
repsPerExperiment = 20;
threshCollision = 0.25; % Anything closer than 25cm crashes
numWindowsPerSecond = 1; % Computed as: expData.paramStruct.spotterCam.fps/expData.paramStruct.frameworkWinSize
maxTimeStepsMOTA = floor(150*numWindowsPerSecond); % Up to 300s

%% Pre-Process real experiments (combine multiple drones into single experiment)
if ~skipPreProcess
	% Ours
	for expInd = 2:3 %1:5
		aux = load([realOutputFolder 'Ours/experiment' num2str(expInd) '/Actuation_oursReal_3N_2x2x1.mat']);
		if ~exist('paramStruct','var')
			paramStruct = aux.paramStruct;
			variableStruct = aux.variableStruct;
		else
			variableStruct = [variableStruct aux.variableStruct(end)];
		end
	end
	preProcessedResultsFileName = strjoin([realOutputFolderAndPrefix '_' paramStruct.typeOfMotion '_' paramStruct.N 'N_' string(paramStruct.roomDimensions).join('x') '.mat'], '');
	save(preProcessedResultsFileName, 'paramStruct','variableStruct', '-v7.3');
	dispImproved(sprintf('Saved pre-processed experiment results as "%s"!\n', preProcessedResultsFileName), 'keepthis');

	% All but ours
	paramStruct = rmfield(paramStruct, {'ourActuationNumLowRiskIterations', 'ourActuationNumBestAssignments', 'ourActuationStopLowRiskCriteria'});
	paramStruct.sigmaNoiseMotion = 0.04;
	experimentAssignments = [1:3; 4:6; 7:9; 10:12; 13:15; 1:5:15; 2:5:15; 3:5:15; 4:5:15; 5:5:15];
	for typeOfMotionCell = {'Random', 'Hovering', 'Landed', 'OneAtATime'}
		typeOfMotion = typeOfMotionCell{:};
		paramStruct.typeOfMotion = [lower(typeOfMotion(1)) typeOfMotion(2:end)];
		data = loadRealExperimentData(struct('datetime',strcat('experiment', strsplit(num2str(1:15),' ')), 'ch',''), ['/Users/carlitos/GoogleDrive/UNI CARLOS/Grad school/Research/Mobisys 18 - Paper/Code/data/Real/' typeOfMotion]);
		
		variableStruct(1:end) = [];
		for experimentInd = 1:size(experimentAssignments,1)
			runningCorrStruct = runMatchingFrameworkOnGivenData(data(experimentAssignments(experimentInd,:)), [],[], paramStruct.frameworkWinSize);
			variableStruct(experimentInd).experimentRanOutOfTime = true;
			variableStruct(experimentInd).groundTruthAssignment = (1:paramStruct.N)';
			variableStruct(experimentInd).posUAVcam = runningCorrStruct.posCam;
			variableStruct(experimentInd).posUAVgt = runningCorrStruct.posCam;
			for f = {'accelUAV', 'accelCam', 't', 'runningWinScore', 'runningLikelihood', 'runningPosterior', 'assignedMatch'}
				variableStruct(experimentInd).(f{:}) = runningCorrStruct.(f{:});
			end
			% plotExperimentResults(runningCorrStruct); % Optional: plot results
		end
		preProcessedResultsFileName = strjoin([realOutputFolderAndPrefix '_' paramStruct.typeOfMotion '_' paramStruct.N 'N_' string(paramStruct.roomDimensions).join('x') '.mat'], '');
		save(preProcessedResultsFileName, 'paramStruct','variableStruct', '-v7.3');
		dispImproved(sprintf('Saved pre-processed experiment results as "%s"!\n', preProcessedResultsFileName), 'keepthis');
	end
end

%% Process all experiment data saved
if ~skipProcess
	for folderAndPrefixCell = {realOutputFolderAndPrefix, simOutputFolderAndPrefix}
		folderAndPrefix = folderAndPrefixCell{:};
		filesInFolder = dir([folderAndPrefix '_*.mat']);
		
		for roomSizeCell = {'4x4x2.mat', '2x2x1.mat', }
			roomSize = roomSizeCell{:};
			processedResultsFileName = [folderAndPrefix '_processed_' roomSize];

			% Initialize experimentsStruct based on the names of the files found
			experimentsStruct = struct();
			for fileNameCell = {filesInFolder.name}
				fileName = fileNameCell{:};
				fileInfo = strsplit(fileName, '_');
				if strcmpi(fileInfo{2}, 'processed'), continue; end % Skip processed results
				if strcmpi(fileInfo{4}, roomSize) % && ~strcmpi(fileInfo{2}, 'landed') % Only process files of the right type of experiment
					N = str2double(fileInfo{3}(1:end-1));
					if ~isfield(experimentsStruct, fileInfo{2}) % Then, arrange the *.mat's according to their motion type (fileInfo{2}, ie. random, hovering, landed...)
						experimentsStruct.(fileInfo{2}).fileNames = {fileName};	     % Initialize the field if it didn't exist
						experimentsStruct.(fileInfo{2}).N = N;
					else
						experimentsStruct.(fileInfo{2}).fileNames{end+1} = fileName; % Otherwise just append the fileName to the end of the list
						experimentsStruct.(fileInfo{2}).N(end+1) = N;
					end
				end
			end

			%% Actually do the processing
			for typeOfMotionCell = fieldnames(experimentsStruct)'
				typeOfMotion = typeOfMotionCell{:};

				% Sort fileNames by N, and save N as a field
				[experimentsStruct.(typeOfMotion).N, sortedInds] = sort(experimentsStruct.(typeOfMotion).N);
				experimentsStruct.(typeOfMotion).N = unique(experimentsStruct.(typeOfMotion).N); % There might be repeated values of N (different normalizations), so keep only unique values
				experimentsStruct.(typeOfMotion).fileNames = experimentsStruct.(typeOfMotion).fileNames(sortedInds); % Sort fileNames by N

				% Initialize remaining fields in the structure
				experimentsStruct.(typeOfMotion).idMOTA = NaN(repsPerExperiment, length(experimentsStruct.(typeOfMotion).N), length(confidenceThresholds), maxTimeStepsMOTA);
				experimentsStruct.(typeOfMotion).idAccuracy = zeros(repsPerExperiment, length(experimentsStruct.(typeOfMotion).N), length(confidenceThresholds));
				experimentsStruct.(typeOfMotion).idAccuracyFirstEver = zeros(size(experimentsStruct.(typeOfMotion).idAccuracy));
				experimentsStruct.(typeOfMotion).idTime = NaN(size(experimentsStruct.(typeOfMotion).idAccuracy));
				experimentsStruct.(typeOfMotion).correctIdTime = NaN(size(experimentsStruct.(typeOfMotion).idMOTA));
				experimentsStruct.(typeOfMotion).correctOverTime = NaN(size(experimentsStruct.(typeOfMotion).idMOTA));
				experimentsStruct.(typeOfMotion).percentCorrectOverTime = NaN(size(experimentsStruct.(typeOfMotion).idMOTA));
				experimentsStruct.(typeOfMotion).numSurvivors = NaN(repsPerExperiment, length(experimentsStruct.(typeOfMotion).N), maxTimeStepsMOTA);
				experimentsStruct.(typeOfMotion).posteriorToConfidenceMapping = zeros(2, length(experimentsStruct.(typeOfMotion).N), length(confidenceThresholds)-1);

				for fileNameInd = 1:numel(experimentsStruct.(typeOfMotion).fileNames)
					fileName = experimentsStruct.(typeOfMotion).fileNames{fileNameInd};
					expData = load(fileName); % Load the experiment mat
					if fileNameInd==1, experimentsStruct.(typeOfMotion).paramStruct = expData.paramStruct; end % Copy static params so they can be accessed later on :)
					N = expData.paramStruct.N; [~,Nind] = find(experimentsStruct.(typeOfMotion).N==N);

					% Traverse every experiment repetition for these variables
					for experimentInd = 1:length(expData.variableStruct)
						% Old experiments need fix: assignedMatch should be initialized randomly (to make accuracy metrics realistic)
						expData.variableStruct(experimentInd).assignedMatch(:,1:expData.paramStruct.frameworkWinSize,:) = repmat(randperm(N)', 1,expData.paramStruct.frameworkWinSize,length(expData.paramStruct.frameworkWinSize));
						expData.variableStruct(experimentInd).runningPosterior(:,:,1:expData.paramStruct.frameworkWinSize) = 1/N;

						% Compute longest time we can evaluate MOTA and survivors for, for this experimentInd
						actualTimeStepsMOTA = min(floor(expData.variableStruct(experimentInd).t(end)*numWindowsPerSecond), maxTimeStepsMOTA);
						actualTimeIndsMOTA = expData.paramStruct.frameworkWinSize*(0:actualTimeStepsMOTA)+1;

						% Obtain the runningPosterior of all assigned drones (ie, at every instant, from the NxM posterior matrix, extract the M values indicated by assignedMatch)
						validPosterior = NaN(size(expData.variableStruct(experimentInd).assignedMatch));
						assignedPosteriorInds = sub2ind(size(expData.variableStruct(experimentInd).runningPosterior), ...
							repmat(1:N, 1,size(expData.variableStruct(experimentInd).assignedMatch, 2)), ...
							reshape(expData.variableStruct(experimentInd).assignedMatch, 1,[]), ...
							reshape(repmat(1:size(expData.variableStruct(experimentInd).assignedMatch,2), size(expData.variableStruct(experimentInd).assignedMatch,1),1), 1,[]));
						validPosterior(~isnan(assignedPosteriorInds)) = expData.variableStruct(experimentInd).runningPosterior(assignedPosteriorInds(~isnan(assignedPosteriorInds))); % Avoid indexing with NaN indices (crashes) -> Initialize validPosterior with NaNs and only overwrite at non-NaN assignedPosteriorInds

						% Compute time & accuracy stats for every possible confidence threshold
						confidentAssignments = expData.variableStruct(experimentInd).assignedMatch; % Make a copy of assignedMatch so we don't modify the original data
						gtAssignment = repmat(expData.variableStruct(experimentInd).groundTruthAssignment, 1,size(confidentAssignments,2));
						for confidenceThreshInd = 1:length(confidenceThresholds)
							confidenceThresh = confidenceThresholds(confidenceThreshInd);
							confidentAssignments(validPosterior<confidenceThresh) = NaN;
							numConfidentAssignmentsOverTime = sum(~isnan(confidentAssignments), 1);
							correctOverTime = sum(confidentAssignments==gtAssignment, 1, 'omitNaN'); correctOverTime(:,1:expData.paramStruct.frameworkWinSize,:) = 1; % First guess is random -> 1/N "success rate"
							percentCorrectOverTime = correctOverTime ./ max(1, numConfidentAssignmentsOverTime); % Make sure we don't divide 0/0. Make divisor >=1 so output would be 0/1=0 in that case
							indFirstAllIDd = find(numConfidentAssignmentsOverTime == N, 1); % First point in time where all assignments made were above confidenceThresh
							indFirstAllIDdCorrectly = find(correctOverTime == N, 1); % First point in time where all assignments made (above confidenceThresh) were correct
							correctFirstIDdIndiv = false(1,N); % For the first guess ever for each drone (row), note whether it was correct or not
							for i=1:N
								indFirstIDindiv = find(~isnan(confidentAssignments(i,:)), 1);
								if ~isempty(indFirstIDindiv), correctFirstIDdIndiv(i) = (confidentAssignments(i,indFirstIDindiv) == gtAssignment(i,1)); end
							end

							% MOTA: # correct guesses / # total guesses made
							experimentsStruct.(typeOfMotion).idMOTA(experimentInd,Nind,confidenceThreshInd,1:actualTimeStepsMOTA+1) = cumsum(correctOverTime(actualTimeIndsMOTA)) ./ max(1, cumsum(numConfidentAssignmentsOverTime(actualTimeIndsMOTA)));
							%experimentsStruct.(typeOfMotion).idMOTA(experimentInd,Nind,confidenceThreshInd,actualTimeStepsMOTA+2:end) = experimentsStruct.(typeOfMotion).idMOTA(experimentInd,Nind,confidenceThreshInd,actualTimeStepsMOTA+1);
							
							% correctOverTime: How many correct guesses per iteration
							experimentsStruct.(typeOfMotion).correctOverTime(experimentInd,Nind,confidenceThreshInd,1:actualTimeStepsMOTA+1) = correctOverTime(actualTimeIndsMOTA);
							experimentsStruct.(typeOfMotion).percentCorrectOverTime(experimentInd,Nind,confidenceThreshInd,1:actualTimeStepsMOTA+1) = percentCorrectOverTime(actualTimeIndsMOTA);

							% idAccuracyFirstEver: % of correct guesses only taking into account the first guess (above confidenceThresh) ever made for each drone (row) individually
							experimentsStruct.(typeOfMotion).idAccuracyFirstEver(experimentInd,Nind,confidenceThreshInd) = sum(correctFirstIDdIndiv)/N;

							% Accuracy: at the first point in time where all assignments were above the confidenceThresh, how many of them were correct
							% idTime: How long it took the system to first have all assigments above the confidenceThresh
							% Note: It could be possible that (eg, for high confidenceThresh values) the system never fully IDs all drones. idTime is initialized with NaNs so just don't write to it in that case
							if ~isempty(indFirstAllIDd)
								experimentsStruct.(typeOfMotion).idAccuracy(experimentInd,Nind,confidenceThreshInd) = correctOverTime(indFirstAllIDd)/N;
								experimentsStruct.(typeOfMotion).idTime(experimentInd,Nind,confidenceThreshInd) = expData.variableStruct(experimentInd).t(indFirstAllIDd);
							end

							% correctIdTime: How long it took the system to first guess all N drones correctly (above confidenceThresh)
							if ~isempty(indFirstAllIDdCorrectly)
								experimentsStruct.(typeOfMotion).correctIdTime(experimentInd,Nind,confidenceThreshInd) = expData.variableStruct(experimentInd).t(indFirstAllIDdCorrectly);
							end
							
							dispImproved(sprintf('\nProcessing experiment results for motion "%s" and roomSize "%s":  N=%2d (%2d/%2d), experiment=%2d/%2d, confidenceThresh=%3d%%...', typeOfMotion, roomSize, N, fileNameInd, numel(experimentsStruct.(typeOfMotion).fileNames), experimentInd, length(expData.variableStruct), round(100*confidenceThresh)));
						end

						% Update the posterior -> confidence mapping
						actualPosterior = expData.variableStruct(experimentInd).runningPosterior(:,:,1:expData.paramStruct.frameworkWinSize:end);
						actualAssignedMatch = expData.variableStruct(experimentInd).assignedMatch(:,1:expData.paramStruct.frameworkWinSize:end);
						assignedPosteriorInds = sub2ind(size(actualPosterior), repmat(1:N, 1,size(actualPosterior, 3)), reshape(actualAssignedMatch, 1,[]), reshape(repmat(1:size(actualPosterior,3), size(actualAssignedMatch,1),1), 1,[]));
						assignedPosteriors = NaN(size(actualAssignedMatch)); assignedPosteriors(~isnan(assignedPosteriorInds)) = actualPosterior(assignedPosteriorInds(~isnan(assignedPosteriorInds)));
						cntTotal = histcounts(assignedPosteriors, confidenceThresholds);
						incorrectAssignments = actualAssignedMatch~=repmat(expData.variableStruct(experimentInd).groundTruthAssignment, 1,size(actualAssignedMatch,2));
						assignedPosteriors(incorrectAssignments) = NaN;
						cntCorrect = histcounts(assignedPosteriors, confidenceThresholds);
						experimentsStruct.(typeOfMotion).posteriorToConfidenceMapping(:,Nind,:) = experimentsStruct.(typeOfMotion).posteriorToConfidenceMapping(:,Nind,:) + permute([cntCorrect; cntTotal], [1 3 2]);

						% Compute the pairwise distance between all drones at all time instants
						distAmongDrones = zeros(N*(N-1)/2, actualTimeIndsMOTA(end));
						survivingDrones = true(N, actualTimeIndsMOTA(end));
						avoidDiagDist = threshCollision.*eye(N); % Use this helper matrix to avoid finding drone crashes of a drone with itself
						prevSurvivors = true(N,1);
						deadDronesW = [];
						for tInd = 1:actualTimeIndsMOTA(end)
							reshapedPos = reshape(expData.variableStruct(experimentInd).posUAVgt(:,tInd,:), [],size(expData.variableStruct(experimentInd).posUAVgt,3));	% Mx3
							distAmongDrones(:,tInd) = pdist(reshapedPos);
							[deadDronesI, deadDronesJ] = find((squareform(distAmongDrones(:,tInd)) + avoidDiagDist) < threshCollision);
							
							if isequal(folderAndPrefix, simOutputFolderAndPrefix)  % Only account for walls in simulation
								distToWalls = [reshapedPos - zeros(size(reshapedPos)), repmat(expData.paramStruct.roomDimensions, size(reshapedPos,1),1) - reshapedPos];
								distToWalls(:,3) = [];	% Ignore floor crashes
								deadDronesW = find(any(distToWalls < 0.7*threshCollision, 2));
							end
							
							shouldDie = unique([deadDronesI; deadDronesJ; deadDronesW]);
							shouldDie(prevSurvivors(shouldDie)==false) = []; % Dead drones can't "kill" other drones
							prevSurvivors(shouldDie) = 0;
							survivingDrones(:,tInd) = prevSurvivors;
						end
						experimentsStruct.(typeOfMotion).numSurvivors(experimentInd,Nind,1:actualTimeStepsMOTA+1) = sum(survivingDrones(:,actualTimeIndsMOTA),1)/N;
	% 					[distHistCount, distHistCenters] = hist(distAmongDrones(:), distHistNbins);
	% 					experimentsStruct.(typeOfMotion).distAmongDrones(experimentInd,Nind,:,:) = [distHistCenters', distHistCount'];
					end
				end
			end
			save(processedResultsFileName, 'experimentsStruct', '-v7.3');
			dispImproved(sprintf('Saved processed experiment results as "%s"!\n', processedResultsFileName), 'keepthis');
		end
	end
	dispImproved(sprintf('Done processing experiment results!\n'), 'keepthis');
end

%%
% plotExperimentResults(struct('experimentMatFileName',[simOutputFolderAndPrefix '_random_5N_norm1_5x5x2.5.mat'], 'experimentInd',1), struct('rawAccel',true, 'scatterCorr',false, 'runningLikelihoodVsWinSize',false, 'runningLikelihoodFull',true));
% generateDronesInRoomVideo([], struct('experimentMatFileName',[simOutputFolderAndPrefix '_random_5N_norm1_5x5x2.5.mat'], 'experimentInd',1), [], [], false);

%% Plot sensor data processing figure
if false
	data = loadRealExperimentData(struct('datetime',{'2017-02-19 17-56-48'}, 'ch','75')); % '2017-02-19 17-59-23','2017-02-19 18-01-29','2017-02-19 18-22-23','2017-02-19 18-24-25'
	d = 3; strAx = char('X'+d-1); % z-axis
	tLims = [0 40]; %data.a_cam.(strAx).t(end) - [35 5];
	tIndsCam = find(tLims(1) <= data.a_cam.(strAx).t & data.a_cam.(strAx).t <= tLims(2));
	tIndsUAV = find(tLims(1) <= data.a_UAV.(strAx).t & data.a_UAV.(strAx).t <= tLims(2));
	tIndsUAVorig = find(tLims(1) <= data.a_UAV_orig.(strAx).t & data.a_UAV_orig.(strAx).t <= tLims(2));
	lWidth = 1.5; fSizeLabels = 15.5; fSizeAxes = 13.5;

	ax = gobjects(2,2);
	figure('Units','pixels', 'Position',[200 200, 1000 230]);
	ax(1,1) = subplot(2,2,1); plot(data.a_UAV_orig.(strAx).t(tIndsUAVorig)-tLims(1), data.a_UAV_orig.(strAx).measured(tIndsUAVorig), 'LineWidth',lWidth);
	ax(1,2) = subplot(2,2,2); plot(data.a_UAV.(strAx).t(tIndsUAV)-tLims(1), data.a_UAV.(strAx).measured(tIndsUAV), 'LineWidth',lWidth);
	ylabel('Accel. (m/s^2)  ', 'FontSize',fSizeLabels);
	ax(2,1) = subplot(2,2,3); plot(data.p_cam.(strAx).t(tIndsCam)-tLims(1), data.p_cam.(strAx).measured(tIndsCam) - min(data.p_cam.(strAx).measured), 'LineWidth',lWidth);
	yticks([0 1]);
	ax(2,2) = subplot(2,2,4); plot(data.a_cam.(strAx).t(tIndsCam)-tLims(1), data.a_cam.(strAx).measured(tIndsCam), 'LineWidth',lWidth);
	ylabel('Accel. (m/s^2)  ', 'FontSize',fSizeLabels);

	xlim(ax(:), [0 diff(tLims)]); ylim([ax(1,:) ax(2,2)], 1.4.*[-1 1]); set(ax(:), 'Box','on', 'FontSize',fSizeAxes);
	for i=1:2, for j=1:2, set(ax(i,j), 'Position', get(ax(i,j),'Position')+[0 0.07-0.01*i 0 -0.1]); end; end
	suplabel('{\bfa_{D,raw}} and {\bfa_{D,filt}}:  Acceleration from worker''s{\bfon-board IMU}', 't', ax(1,:), 0, 'FontWeight','normal', 'FontSize',fSizeLabels+1.5);
	suplabel('{\bfx_E\rightarrowa_E}:  Position\rightarrowAcceleration from spotter''s{\bfcamera}', 't', ax(2,:), 0, 'FontWeight','normal', 'FontSize',fSizeLabels+1.5);
	suplabel('Accel. (m/s^2)  ', 'y', ax(1,:), 0, 'FontSize',fSizeLabels); suplabel('Pos. (m)', 'y', ax(2,:), 0, 'FontSize',fSizeLabels); suplabel('Time (s)', 'x', ax(:,1), 0, 'FontSize',fSizeLabels); suplabel('Time (s)', 'x', ax(:,2), 0, 'FontSize',fSizeLabels);
	saveFigToFile('aCam_vs_aUAV');
end

%% Plot iterative update figure
if false
	data = loadRealExperimentData(struct('datetime',{'2017-02-19 17-59-23','2017-02-19 18-01-29','2017-02-19 18-22-23','2017-02-19 18-24-25'}, 'ch','75')); % '2017-02-19 17-59-23','2017-02-19 18-01-29','2017-02-19 18-22-23','2017-02-19 18-24-25'
	d = 3; strAx = char('X'+d-1); % z-axis
	runningCorrStruct = runMatchingFrameworkOnGivenData(data, [], [], 12, d);
	tLims = [0 30]; %data.a_cam.(strAx).t(end) - [35 5];
	tIndsCam = cell(1,length(data)); tIndsUAV = cell(1,length(data));
	for i = 1:length(data)
		tIndsCam{i} = find(tLims(1) <= data(i).a_cam.(strAx).t & data(i).a_cam.(strAx).t <= tLims(2));
		tIndsUAV{i} = find(tLims(1) <= data(i).a_UAV.(strAx).t & data(i).a_UAV.(strAx).t <= tLims(2));
	end
	lWidth = 2.5; fSizeLabels = 15; fSizeAxes = 13.5;

	ax = gobjects(2,length(data));
	figure('Units','pixels', 'Position',[200 200, 1000 210]);
	iCam = 1;
	for i = 1:length(data)
		ax(1,i) = subplot(2,length(data),i); plot(data(i).a_UAV.(strAx).t(tIndsUAV{i})-tLims(1), data(i).a_UAV.(strAx).measured(tIndsUAV{i}), 'LineWidth',lWidth);
		hold on; plot(data(iCam).a_cam.(strAx).t(tIndsCam{iCam})-tLims(1), data(iCam).a_cam.(strAx).measured(tIndsCam{iCam}), ':', 'Color',[0.83,0,0.1], 'LineWidth',lWidth*1.25);

		ax(2,i) = subplot(2,length(data),i+length(data)); plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningPosterior(i,iCam,:)), 'LineWidth',2, 'Color',[0.93,0.69,0.125]);
	end

	xlim(ax(:), [0 diff(tLims)]); ylim(ax(1,:), 1.*[-1 1]); ylim(ax(2,:), [0 1]); yticks(ax(2,:), [0,1]); set(ax(:), 'Box','on', 'FontSize',fSizeAxes);
	for i=1:2, for j=1:length(data), set(ax(i,j), 'Position', get(ax(i,j),'Position')+[0 0.08 0 -0.02]); end; end
	for i =1:length(data), suplabel('Time (s)', 'x', ax(:,i), 0, 'FontSize',fSizeLabels); end
	suplabel({'Accel.','(m/s^2)'}, 'y', ax(1,:), 0.005, 'FontSize',fSizeLabels); suplabel({'Posterior','P(\theta_{i,j}|D^t)'}, 'y', ax(2,:), 0.0, 'FontSize',fSizeLabels);
	l=legend(ax(1,end), {' {\bfa}_{Dj}', ' {\bfa}_{E1}'}, 'Orientation','vertical', 'FontSize',fSizeAxes-1); set(l, 'Units','normalized', 'Position',[0.926 0.8 0.045 0.05]);
	saveFigToFile('iterativeUpdate');
end



%% Plot accuracy vs TYPE OF MOTION over time (evaluation of Matching only - no actuation yet)
if true
	tEnd = 20; intervalErrorBar = 1;
	lWidth = 2; lWidthErr = 1; fSizeLabels = 17; fSizeAxes = 13.5;
	typeOfMotionCell = {'random', 'hovering', 'landed'};
	legendCell = {'Random motion', 'Hovering', 'Landed'};
	for i = 1:3
		if i <= 2, N = 3; roomSize = '2x2x1.mat'; else, N = 24; roomSize = '4x4x2.mat'; end
		if i <= 1, folderAndPrefix = realOutputFolderAndPrefix; else, folderAndPrefix = simOutputFolderAndPrefix; end
		if isequal(folderAndPrefix, simOutputFolderAndPrefix), simOrReal = 'Sim'; else, simOrReal = 'Real'; end
		processedResultsFileName = [folderAndPrefix '_processed_' roomSize];
		load(processedResultsFileName);

		figure('Units','pixels', 'Position',[200 200, 560 200]); hold on;
		h = gobjects(1,numel(typeOfMotionCell));
		for typeOfMotionInd = 1:numel(typeOfMotionCell)
			typeOfMotion = typeOfMotionCell{typeOfMotionInd};
			[~,Nind] = find(experimentsStruct.(typeOfMotion).N==N);

			% Compute accuracy and numSurvivors stats (mean, std)
			[MOTAmean, MOTAstd, numSurvivorsMean, numSurvivorsStd] = computeMotaAndSurvivalStats(experimentsStruct, typeOfMotion, intervalErrorBar);
			t = (0:size(numSurvivorsMean,2)-1) / numWindowsPerSecond;

			h(typeOfMotionInd) = errorbar(t, MOTAmean(Nind,:), MOTAstd(Nind,:), 'LineWidth',lWidthErr);
			plot(t, MOTAmean(Nind,:), 'LineWidth',lWidth, 'Color',get(h(typeOfMotionInd),'Color'));
		end
		xlim([0,tEnd]); ylim([0 100]);
		xlabel('Time (s)', 'FontSize',fSizeLabels); ylabel('Accuracy (%)', 'FontSize',fSizeLabels+1);
		set(gca, 'Box','on', 'FontSize',fSizeAxes);
		l=legend(h, legendCell, 'FontSize',fSizeAxes, 'Location','NorthEast');% set(l, 'Units','normalized', 'Position',[0.17 0.494 0.7 0]);
		saveFigToFile(['AccuracyOverTimeForDifferentTypeOfMotion_' num2str(N) 'N_' simOrReal]);
	end
end

%% Plot accuracy and survival rate vs TYPE OF MOTION over time (actuation)
if true
	tEnd = 20; intervalErrorBar = 2;
	for i = 1:3
		if i <= 2, N = 3; roomSize = '2x2x1.mat'; else, N = 24; roomSize = '4x4x2.mat'; end
		if i <= 1, folderAndPrefix = realOutputFolderAndPrefix; else, folderAndPrefix = simOutputFolderAndPrefix; end
		if isequal(folderAndPrefix, simOutputFolderAndPrefix), simOrReal = 'Sim'; else, simOrReal = 'Real'; end

		processedResultsFileName = [folderAndPrefix '_processed_' roomSize];
		load(processedResultsFileName);
		typeOfMotionCell = {'random', 'oneAtATime', 'ours'};
		legendCell = {'Random', 'One-by-one', '\textit{IDrone} (ours)'};
		typeOfMotionCell{3} = [typeOfMotionCell{3} simOrReal];

		plotAccuracyAndSurvivorVsTypeOfMotion(experimentsStruct, typeOfMotionCell, legendCell, N, numWindowsPerSecond, tEnd, intervalErrorBar);
		saveFigToFile(['AccuracyAndSurvivalOverTimeForDifferentTypeOfActuation_' num2str(N) 'N_' simOrReal]);
	end
end

%% Plot MOTA and survivors vs DRONE DENSITY over time
if true
	tEnd = 30; intervalErrorBar = 2;
	for roomSizeCell = {'4x4x2.mat'}
		roomSize = roomSizeCell{:};
		processedResultsFileName = [simOutputFolderAndPrefix '_processed_' roomSize];
		load(processedResultsFileName);
		typeOfMotionCell = fieldnames(experimentsStruct)';

		for typeOfMotionInd = 1:numel(typeOfMotionCell)
			typeOfMotion = typeOfMotionCell{typeOfMotionInd};

			plotAccuracyAndSurvivorVsN(experimentsStruct, typeOfMotion, numWindowsPerSecond, tEnd, intervalErrorBar);
			saveFigToFile(['AccuracyAndSurvivalOverTimeForDifferentN_' typeOfMotion]);
		end
	end
end



%% Plot posterior vs confidence
if false
	lWidth = 2; fSizeLabels = 17; fSizeAxes = 14;
	tEnd = 100;
	for roomSizeCell = {'5x5x2.5.mat'} %{'5x5x2.5.mat', '15x10x3.mat'}
		roomSize = roomSizeCell{:};
		processedResultsFileName = [simOutputFolderAndPrefix '_processed_' roomSize];
		load(processedResultsFileName);
		typeOfMotionCell = fieldnames(experimentsStruct)';
		N = 5:5:20;
		cntCorrect = zeros(length(N), length(confidenceThresholds)-1);
		cntTotal = zeros(size(cntCorrect));

		% Compute accuracy
		for NarrInd = 1:numel(N)
			Nind = N(NarrInd);%-1;
			for typeOfMotionInd = 1:numel(typeOfMotionCell)
				typeOfMotion = typeOfMotionCell{typeOfMotionInd};
				cntCorrect(NarrInd,:) = cntCorrect(NarrInd,:) + reshape(experimentsStruct.(typeOfMotion).posteriorToConfidenceMapping(1,Nind,:), 1,[]);
				cntTotal(NarrInd,:) = cntTotal(NarrInd,:) + reshape(experimentsStruct.(typeOfMotion).posteriorToConfidenceMapping(2,Nind,:), 1,[]);
			end
		end
		cntTotal(cntTotal<30)=0;
		
		% Plot results
		figure('Units','pixels', 'Position',[200 200, 560 140]); hold on;
		plot(confidenceThresholds(1:end-1)+0.5*diff(confidenceThresholds), cntCorrect./cntTotal, 'LineWidth',lWidth);
		legend(cellstr(strcat('N=',num2str((N)'))), 'Orientation','vertical', 'FontSize',fSizeAxes, 'Location','SouthEast');

		xlim([0,1]); ylim([0 1]); xlabel('Posterior (%)', 'FontSize',fSizeLabels+1); ylabel('Confidence (%)', 'FontSize',fSizeLabels+1);
		set(gca, 'Box','on', 'FontSize',fSizeAxes);
		saveFigToFile(['posteriorVsConfidence']);
	end
end

%% Plot MOTA and survivors vs different actuation switches over time
if false
	lWidth = 3; lWidthErr = 1; fSizeLabels = 17; fSizeAxes = 13.5;
	tEnd = 100;
	N = 20;
	Q = [2 4 8 16 32];
	roomSize = '5x5x2.5.mat';
	typeOfMotion = 'hovering'; %'ours';
	processedResultsFileName = [simOutputFolderAndPrefix '_processed_' roomSize];
	load(processedResultsFileName);
	
	figure('Units','pixels', 'Position',[200 200, 560 375]);
	ax = gobjects(2,1); h = gobjects(2,numel(Q));

	% Compute MOTA and numSurvivors stats (mean, std)
	MOTAmean = squeeze(mean(experimentsStruct.(typeOfMotion).idMOTA(:,:,1,:), 1, 'omitnan')); MOTAmean(:)=0;
	MOTAstd = squeeze(std(experimentsStruct.(typeOfMotion).idMOTA(:,:,1,:), 0,1, 'omitnan')); MOTAstd(:)=0;
	numSurvivorsMean = squeeze(mean(experimentsStruct.(typeOfMotion).numSurvivors, 1, 'omitnan')); numSurvivorsMean(:)=0;
	numSurvivorsStd = squeeze(std(experimentsStruct.(typeOfMotion).numSurvivors, 0,1, 'omitnan')); numSurvivorsStd(:)=0;
	t = (0:size(numSurvivorsMean,2)-1) / numWindowsPerSecond;

	for Qind = 1:numel(Q)

		% Plot MOTA
		ax(1) = subplot(2,1,1); hold on;
% 		errorbar(t, MOTAmean(Nind,:), MOTAstd(Nind,:), 'LineWidth',lWidthErr);
		h(1,Qind) = plot(t, MOTAmean(Qind,:), 'LineWidth',lWidth); %, 'Color',get(h(1,NarrInd),'Color'));
		xlim([0,tEnd]); ylim([0 1]); ylabel('MOTA (%)', 'FontSize',fSizeLabels+1);

		% Plot survival rate
		ax(2) = subplot(2,1,2); hold on;
% 		h(2,NarrInd) = errorbar(t, numSurvivorsMean(Nind,:), numSurvivorsStd(Nind,:), 'LineWidth',lWidthErr);
% 		plot(t, numSurvivorsMean(Nind,:), 'LineWidth',lWidth, 'Color',get(h(2,NarrInd),'Color'));
		xlim([0,tEnd]); ylim([-0.05 1.05]); ylabel('Survival rate (%)', 'FontSize',fSizeLabels+1);
	end

	suplabel('Time (s)', 'x', ax(:), -0.01, 'FontSize',fSizeLabels);
	set(ax(:), 'Box','on', 'FontSize',fSizeAxes);
 	l=legend(h(1,:), cellstr(strcat('Q=',num2str((Q)'))), 'Orientation','horizontal', 'FontSize',fSizeAxes); set(l, 'Units','normalized', 'Position',[0.17 0.494 0.7 0]);
	saveFigToFile('MOTAandSurvivalOverTimeForDifferentQ');
end

%% Plot raw accel for real vs simulation
if false
	data = loadRealExperimentData(struct('datetime',{'2017-02-19 17-59-23','2017-02-19 18-01-29','2017-02-19 18-22-23','2017-02-19 18-24-25'}, 'ch','75')); % '2017-02-19 17-59-23','2017-02-19 18-01-29','2017-02-19 18-22-23','2017-02-19 18-24-25'
	d = 3; strAx = char('X'+d-1); % z-axis
	runningCorrStruct = runMatchingFrameworkOnGivenData(data, [], [], 12, d);
	tLims = [0 30]; %data.a_cam.(strAx).t(end) - [35 5];
	tIndsCam = cell(1,length(data)); tIndsUAV = cell(1,length(data));
	for i = 1:length(data)
		tIndsCam{i} = find(tLims(1) <= data(i).a_cam.(strAx).t & data(i).a_cam.(strAx).t <= tLims(2));
		tIndsUAV{i} = find(tLims(1) <= data(i).a_UAV.(strAx).t & data(i).a_UAV.(strAx).t <= tLims(2));
	end
	N = 4; Nind = N-1;
	roomSize = '5x5x2.5.mat';
	typeOfMotion = 'random';
	aux = load([simOutputFolderAndPrefix '_simulation_4N_norm0_8x5x2.5.mat']);
	lWidth = 2.5; fSizeLabels = 15; fSizeAxes = 13.5;

	ax = gobjects(2,length(data));
	figure('Units','pixels', 'Position',[200 200, 1000 210]);
	iCam = 1;
	for i = 1:length(data)
		ax(1,i) = subplot(2,length(data),i); plot(data(i).a_UAV.(strAx).t(tIndsUAV{i})-tLims(1), data(i).a_UAV.(strAx).measured(tIndsUAV{i}), 'LineWidth',lWidth);
		hold on; plot(data(i).a_cam.(strAx).t(tIndsCam{i})-tLims(1), data(i).a_cam.(strAx).measured(tIndsCam{i}), 'Color',[0.83,0,0.1], 'LineWidth',lWidth*0.8);

		ax(2,i) = subplot(2,length(data),i+length(data)); plot(aux.variableStruct(1).t, aux.variableStruct(1).accelUAV(i,:,d), 'LineWidth',lWidth);
		hold on; plot(aux.variableStruct(1).t, aux.variableStruct(1).accelCam(i,:,d), '-', 'Color',[0.83,0,0.1], 'LineWidth',lWidth*0.8);
	end

	xlim(ax(:), [0 diff(tLims)]); ylim(ax(:), 1.*[-1 1]); set(ax(:), 'Box','on', 'FontSize',fSizeAxes);
	for i=1:2, for j=1:length(data), set(ax(i,j), 'Position', get(ax(i,j),'Position')+[0 0.08 0 -0.02]); end; end
	for i =1:length(data), suplabel('Time (s)', 'x', ax(:,i), 0, 'FontSize',fSizeLabels); end
	suplabel({'REAL','accel. (m/s^2)'}, 'y', ax(1,:), 0.005, 'FontSize',fSizeLabels); suplabel({'SIMULATION','accel. (m/s^2)'}, 'y', ax(2,:), 0.005, 'FontSize',fSizeLabels);
	l=legend(ax(1,end), {' {\bfa}_{Dj}', ' {\bfa}_{E1}'}, 'Orientation','vertical', 'FontSize',fSizeAxes-1); set(l, 'Units','normalized', 'Position',[0.926 0.8 0.045 0.05]);
	l=legend(ax(2,end), {' {\bfa}_{Dj}', ' {\bfa}_{E1}'}, 'Orientation','vertical', 'FontSize',fSizeAxes-1); set(l, 'Units','normalized', 'Position',[0.926 0.32 0.045 0.05]);
	saveFigToFile('rawAccelRealVsSim');
end

%% Plot MOTA over time for real vs simulation
if false
	lWidth = 3; lWidthErr = 1; fSizeLabels = 17; fSizeAxes = 13.5;
	tEnd = 100;
	N = 5; Nind = N-1;
	roomSize = '5x5x2.5.mat';
	typeOfMotion = 'random';
	processedResultsFileName = [simOutputFolderAndPrefix '_processed_' roomSize];
	load(processedResultsFileName);
	
	figure('Units','pixels', 'Position',[200 200, 560 200]);

	% Compute MOTA and numSurvivors stats (mean, std)
	MOTAmean = squeeze(mean(experimentsStruct.(typeOfMotion).idMOTA(:,:,1,:), 1, 'omitnan')); MOTAmean(:)=0;
	MOTAstd = squeeze(std(experimentsStruct.(typeOfMotion).idMOTA(:,:,1,:), 0,1, 'omitnan')); MOTAstd(:)=0;
	t = (0:size(MOTAmean,4)-1) / numWindowsPerSecond;

	hold on;
% 	errorbar(t, MOTAmean(Nind,:), MOTAstd(Nind,:), 'LineWidth',lWidthErr);
	plot(t, MOTAmean(Nind,:), 'LineWidth',lWidth); %, 'Color',get(h(1,NarrInd),'Color'));
	xlim([0,tEnd]); ylim([0 1]); ylabel('MOTA (%)', 'FontSize',fSizeLabels+1);

	xlabel('Time (s)', 'FontSize',fSizeLabels);
	set(gca, 'Box','on', 'FontSize',fSizeAxes);
	legend({'Real experiments', 'Simulation'}, 'FontSize',fSizeAxes, 'Location','SouthEast');
%  	l=legend(h(1,:), cellstr(strcat('Q=',num2str((Q)'))), 'Orientation','horizontal', 'FontSize',fSizeAxes); set(l, 'Units','normalized', 'Position',[0.17 0.494 0.7 0]);
	saveFigToFile('MOTAForRealVsSim');
end
return

%% Plot figures from processed *.mat files
for roomSizeCell = {'4x4x2.mat'} %{'5x5x2.5.mat', '15x10x3.mat'}
	roomSize = roomSizeCell{:};
	processedResultsFileName = [simOutputFolderAndPrefix '_processed_' roomSize];
	load(processedResultsFileName);
	%idMOTA, idAccuracyFirstEver, [idAccuracy], idTime, [correctIdTime], distAmongDrones
	%%
	% yyaxis left;
	% boxplot(squeeze(experimentsStruct.(typeOfMotion).idTime(:,end,1,:)), confidenceThresholds, 'boxstyle','filled', 'medianstyle','target', 'sym','r.', 'positions',confidenceThresholds);
	% yyaxis right;
	% boxplot(squeeze(experimentsStruct.(typeOfMotion).idMOTA(:,end,1,:)), confidenceThresholds, 'boxstyle','filled', 'medianstyle','target', 'sym','r.', 'positions',confidenceThresholds);
	for N = ceil(24./[5:-1:1])
		Nind = N-1;
		for normTypeInd = 1 %:normalizationTypes
			for typeOfMotionCell = {'random', 'hovering'} % fieldnames(experimentsStruct)'
				typeOfMotion = typeOfMotionCell{:};
				if strcmpi(typeOfMotion, 'random'), tit = 'Random motion'; else, tit = typeOfMotion; end
				meanAccuracy = mean(squeeze(experimentsStruct.(typeOfMotion).idAccuracyFirstEver(:,Nind,normTypeInd,:)), 1);
				meanTime = mean(squeeze(experimentsStruct.(typeOfMotion).idTime(:,Nind,normTypeInd,:)), 1);
				meanCorrectTime = mean(squeeze(experimentsStruct.(typeOfMotion).correctIdTime(:,Nind,normTypeInd,:)), 1);
				ind100percentCorrect = find(meanAccuracy==1, 1);
				figure;
				yyaxis left;
				plot(100*confidenceThresholds(1:end-1), meanAccuracy(1:end-1), 'LineWidth',3); hold on;
				title([tit ', N=' num2str(N) ', norm=' num2str(normTypeInd-1)], 'FontSize',18); xlabel('Confidence threshold (%)', 'FontSize',18); ylabel('ID accuracy (%)', 'FontSize',18); ylim([0, 1]); set(gca, 'FontSize',14);
				% plot([100/N 100/N], [0 1], '--k'); annotation('textarrow', axesCoordsToNormFigCoords('x', 100/N) + [0.1 0], axesCoordsToNormFigCoords('y', 0.5) - [0.1 0], 'String','1/N', 'FontSize',16);
				if ~isempty(ind100percentCorrect)
					halfDx=0.015/2; annotation('rectangle', axesCoordsToNormFigCoords([100*(confidenceThresholds(ind100percentCorrect)-halfDx), 0, 200*halfDx, 1]), 'EdgeColor','green', 'FaceColor','green', 'FaceAlpha',0.2);
				end
				yyaxis right;
				plot(100*confidenceThresholds, meanTime, 'LineWidth',3); hold on;
				plot(100*confidenceThresholds, meanCorrectTime, '--g', 'LineWidth',3);
				ylabel('ID time (s)', 'FontSize',18); set(gca, 'FontSize',14);
				if ~isempty(ind100percentCorrect)
					annotation('textarrow', axesCoordsToNormFigCoords('x', 100*(confidenceThresholds(ind100percentCorrect))) + [0.1 0], ...
						axesCoordsToNormFigCoords('y', meanTime(ind100percentCorrect)) - [0 0], ...
						'String',sprintf('100%% accuracy at %d%% (t=%.2fs)', round(100*confidenceThresholds(ind100percentCorrect)), meanTime(ind100percentCorrect)), 'FontSize',16);
				end
				saveas(gcf, ['idTimeAndAccuracyVsConfidenceThresh_' typeOfMotion '_N' num2str(N) 'norm' num2str(normTypeInd-1)], 'epsc');
			end
		end
	end
	
	%%
	normTypeInd = 1;
	confidenceThreshInd = 1;
	for yAxisCell = {'idAccuracy', 'idMOTA', 'idAccuracyFirstEver', 'idAccuracy', 'idTime', 'correctIdTime'}
		yAxis = yAxisCell{:};
		figure;
		boxplot(experimentsStruct.random.(yAxis)(:,:,normTypeInd,confidenceThreshInd), experimentsStruct.random.N, 'Colors','b');
		hold on;
		boxplot(experimentsStruct.hovering.(yAxis)(:,:,normTypeInd,confidenceThreshInd), experimentsStruct.hovering.N, 'Colors','r');
		title(yAxis); ylabel('ID accuracy (%)', 'FontSize',18); set(gca, 'FontSize',14); %ylim([0, 1]);
		% saveas(gcf, ['idMOTAvsTypeOfMotion_N' num2str(N) excludeLanded], 'epsc');
	end

	%%
	normTypeInd = 1;
	for N = 5:5:25
		Nind = N-1;
		for excludeLandedCell = {'', '_noLanded'}
			excludeLanded = excludeLandedCell{:};
			
			confidenceThreshInd = 16;
			yAxis = 'idAccuracy';
			if isempty(excludeLanded)
				if N>15 || strcmpi(roomSize, '15x10x3.mat'), continue; end % No 'landed' data for N>15
				figure;
				boxplot([experimentsStruct.random.(yAxis)(:,Nind,normTypeInd,confidenceThreshInd), experimentsStruct.hovering.(yAxis)(:,Nind,normTypeInd,confidenceThreshInd), experimentsStruct.landed.(yAxis)(:,Nind,normTypeInd,confidenceThreshInd)], {'Random motion', 'Hovering', 'Landed'});
				hold on; plot([-1 4], repmat(1/N, 1,2), ':k');
			else
				figure;
				boxplot([experimentsStruct.random.(yAxis)(:,Nind,normTypeInd,confidenceThreshInd), experimentsStruct.hovering.(yAxis)(:,Nind,normTypeInd,confidenceThreshInd)], {'Random motion', 'Hovering'});
			end
			ylabel('ID accuracy (%)', 'FontSize',18); title(['N=' num2str(N)], 'FontSize',18); set(gca, 'FontSize',14); %ylim([0, 1]);
			saveas(gcf, ['idMOTAvsTypeOfMotion_N' num2str(N) excludeLanded], 'epsc');

			figure;
			yAxis = 'correctIdTime';
			if isempty(excludeLanded)
				boxplot([experimentsStruct.random.(yAxis)(:,Nind,normTypeInd,confidenceThreshInd), experimentsStruct.hovering.(yAxis)(:,Nind,normTypeInd,confidenceThreshInd), experimentsStruct.landed.(yAxis)(:,Nind,normTypeInd,confidenceThreshInd)], {'Random motion', 'Hovering', 'Landed'});
			else
				boxplot([experimentsStruct.random.(yAxis)(:,Nind,normTypeInd,confidenceThreshInd), experimentsStruct.hovering.(yAxis)(:,Nind,normTypeInd,confidenceThreshInd)], {'Random motion', 'Hovering'});
			end
			ylabel('ID time (s)', 'FontSize',18); title(['N=' num2str(N)], 'FontSize',18); set(gca, 'FontSize',14);
			saveas(gcf, ['idTimeVsTypeOfMotion_N' num2str(N) excludeLanded], 'epsc');
		end
	end

	%%
	for N = [5 10 15]
		Nind = N-1;
		for typeOfMotionCell = fieldnames(experimentsStruct)'
			typeOfMotion = typeOfMotionCell{:};
			if strcmpi(typeOfMotion, 'random'), tit = 'Random motion'; else, tit = typeOfMotion; end
			figure;
			bar(squeeze(experimentsStruct.(typeOfMotion).distAmongDrones(1,Nind,normTypeInd,:,1)), squeeze(experimentsStruct.(typeOfMotion).distAmongDrones(1,Nind,normTypeInd,:,2)));
			title([tit ', N=' num2str(N)], 'FontSize',18); xlabel('Distance among drones (m)', 'FontSize',18); ylabel('Count (#)', 'FontSize',18); set(gca, 'FontSize',14);
			saveas(gcf, ['histDistAmongDrones_' typeOfMotion '_N' num2str(N)], 'epsc');
		end
	end
end