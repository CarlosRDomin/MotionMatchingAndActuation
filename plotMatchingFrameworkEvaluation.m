%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots all figures related to the matching framework evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all;
cdToThisScriptsDirectory;

typeOfExperiment = 'MatchingFramework';
simOutputFolder = [fileparts(mfilename('fullpath')) '/data/Simulations/'];
skipProcess = true; % Whether or not to (re)process all the experiment data

% Constants
confidenceThresholds = 0:0.01:1;
repsPerExperiment = 20;
normalizationTypes = 3;
distHistNbins = 100;

%% Process all experiment data saved
filesInFolder = dir([simOutputFolder typeOfExperiment '_*.mat']);
if ~skipProcess
	for roomSizeCell = {'5x5x2.5.mat', '15x10x3.mat'}
		roomSize = roomSizeCell{:};
		processedResultsFileName = [simOutputFolder typeOfExperiment '_processed_' roomSize];

		% Initialize experimentsStruct based on the names of the files found
		experimentsStruct = struct();
		for fileNameCell = {filesInFolder.name}
			fileName = fileNameCell{:};
			fileInfo = strsplit(fileName, '_');
			if strcmp(fileInfo{2}, 'processed'), continue; end % Skip processed results
			if strcmp(fileInfo{1}, typeOfExperiment) && strcmp(fileInfo{5}, roomSize) % Only process files of the right type of experiment
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

			% Sort fileNames by N, and save N as a 
			[experimentsStruct.(typeOfMotion).N, sortedInds] = sort(experimentsStruct.(typeOfMotion).N);
			experimentsStruct.(typeOfMotion).N = unique(experimentsStruct.(typeOfMotion).N); % There might be repeated values of N (different normalizations), so keep only unique values
			experimentsStruct.(typeOfMotion).fileNames = experimentsStruct.(typeOfMotion).fileNames(sortedInds); % Sort fileNames by N

			% Initialize remaining fields in the structure
			experimentsStruct.(typeOfMotion).idMOTA = NaN(repsPerExperiment, length(experimentsStruct.(typeOfMotion).N), normalizationTypes, length(confidenceThresholds));
			experimentsStruct.(typeOfMotion).idAccuracyFirstEver = zeros(size(experimentsStruct.(typeOfMotion).idMOTA));
			experimentsStruct.(typeOfMotion).idAccuracy = zeros(size(experimentsStruct.(typeOfMotion).idMOTA));
			experimentsStruct.(typeOfMotion).idTime = NaN(size(experimentsStruct.(typeOfMotion).idMOTA));
			experimentsStruct.(typeOfMotion).correctIdTime = NaN(size(experimentsStruct.(typeOfMotion).idMOTA));
			experimentsStruct.(typeOfMotion).distAmongDrones = zeros(repsPerExperiment, length(experimentsStruct.(typeOfMotion).N), normalizationTypes, distHistNbins, 2);

			for fileNameInd = 1:numel(experimentsStruct.(typeOfMotion).fileNames)
				fileName = experimentsStruct.(typeOfMotion).fileNames{fileNameInd};
				expData = load(fileName); % Load the experiment mat
				N = expData.paramStruct.N; [~,Nind] = find(experimentsStruct.(typeOfMotion).N==N);
				fileInfo = strsplit(fileName, '_'); normType = str2double(fileInfo{4}(end)); normTypeInd = normType+1; % Forgot to save normType into the experiment parameters, so extract it from file name ;)

				% Traverse every experiment repetition for these variables
				for experimentInd = 1:length(expData.variableStruct)
					% Old experiments need fix: assignedMatch should be initialized randomly (to make accuracy metrics realistic)
					expData.variableStruct(experimentInd).assignedMatch(:,1:expData.paramStruct.frameworkWinSize,:) = repmat(randperm(N)', 1,expData.paramStruct.frameworkWinSize,length(expData.paramStruct.frameworkWinSize));
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
						correctOverTime = sum(confidentAssignments==gtAssignment, 1, 'omitNaN');
						percentCorrectOverTime = correctOverTime ./ max(1, numConfidentAssignmentsOverTime); % Make sure we don't divide 0/0. Make divisor >=1 so output would be 0/1=0 in that case
						indFirstAllIDd = find(numConfidentAssignmentsOverTime == N, 1); % First point in time where all assignments made were above confidenceThresh
						indFirstAllIDdCorrectly = find(correctOverTime == N, 1); % First point in time where all assignments made (above confidenceThresh) were correct
						correctFirstIDdIndiv = false(1,N); % For the first guess ever for each drone (row), note whether it was correct or not
						for i=1:N
							indFirstIDindiv = find(~isnan(confidentAssignments(i,:)), 1);
							if ~isempty(indFirstIDindiv), correctFirstIDdIndiv(i) = (confidentAssignments(i,indFirstIDindiv) == gtAssignment(i,1)); end
						end

						% MOTA: # correct guesses / # total guesses made
						experimentsStruct.(typeOfMotion).idMOTA(experimentInd,Nind,normTypeInd,confidenceThreshInd) = sum(correctOverTime(:)) / max(1, sum(numConfidentAssignmentsOverTime(:)));

						% idAccuracyFirstEver: % of correct guesses only taking into account the first guess (above confidenceThresh) ever made for each drone (row) individually
						experimentsStruct.(typeOfMotion).idAccuracyFirstEver(experimentInd,Nind,normTypeInd,confidenceThreshInd) = sum(correctFirstIDdIndiv)/N;

						% Accuracy: at the first point in time where all assignments were above the confidenceThresh, how many of them were correct
						% idTime: How long it took the system to first have all assigments above the confidenceThresh
						% Note: It could be possible that (eg, for high confidenceThresh values) the system never fully IDs all drones. idTime is initialized with NaNs so just don't write to it in that case
						if ~isempty(indFirstAllIDd)
							experimentsStruct.(typeOfMotion).idAccuracy(experimentInd,Nind,normTypeInd,confidenceThreshInd) = correctOverTime(indFirstAllIDd)/N;
							experimentsStruct.(typeOfMotion).idTime(experimentInd,Nind,normTypeInd,confidenceThreshInd) = expData.variableStruct(experimentInd).t(indFirstAllIDd);
						end

						% correctIdTime: How long it took the system to first guess all N drones correctly (above confidenceThresh)
						if ~isempty(indFirstAllIDdCorrectly)
							experimentsStruct.(typeOfMotion).correctIdTime(experimentInd,Nind,normTypeInd,confidenceThreshInd) = expData.variableStruct(experimentInd).t(indFirstAllIDdCorrectly);
						end

						dispImproved(sprintf('\nProcessing experiment results for motion "%s" and roomSize "%s":  N=%2d, norm=%d (%2d/%2d), experiment=%2d/%2d, confidenceThresh=%3d%%...', typeOfMotion, roomSize, N, normType, fileNameInd, numel(experimentsStruct.(typeOfMotion).fileNames), experimentInd, length(expData.variableStruct), round(100*confidenceThresh)));
					end

					% Compute the pairwise distance between all drones at all time instants
					distAmongDrones = zeros(N*(N-1)/2, length(expData.variableStruct(experimentInd).t));
					for tInd = 1:length(expData.variableStruct(experimentInd).t)
						distAmongDrones(:,tInd) = pdist(reshape(expData.variableStruct(experimentInd).posUAVgt(:,tInd,:), [],size(expData.variableStruct(experimentInd).posUAVgt,3)));
					end
					[distHistCount, distHistCenters] = hist(distAmongDrones(:), distHistNbins);
					experimentsStruct.(typeOfMotion).distAmongDrones(experimentInd,Nind,normTypeInd,:,:) = [distHistCenters', distHistCount'];
				end
			end
		end
		dispImproved(sprintf('Done processing experiment results!\n'), 'keepthis');
		save(processedResultsFileName, 'experimentsStruct', '-v7.3');
		dispImproved(sprintf('Saved processed experiment results as "%s"!\n', processedResultsFileName), 'keepthis');
	end
end

%%
% plotExperimentResults(struct('experimentMatFileName',[simOutputFolder 'MatchingFramework_random_5N_norm1_5x5x2.5.mat'], 'experimentInd',1), struct('rawAccel',true, 'scatterCorr',false, 'runningLikelihoodVsWinSize',false, 'runningLikelihoodFull',true));
% generateDronesInRoomVideo([], struct('experimentMatFileName',[simOutputFolder 'MatchingFramework_random_5N_norm1_5x5x2.5.mat'], 'experimentInd',1), [], [], false);

%% Plot figures from processed *.mat files
for roomSizeCell = {'5x5x2.5.mat'}%{'5x5x2.5.mat', '15x10x3.mat'}
	roomSize = roomSizeCell{:};
	processedResultsFileName = [simOutputFolder typeOfExperiment '_processed_' roomSize];
	load(processedResultsFileName);
	%idMOTA, idAccuracyFirstEver, [idAccuracy], idTime, [correctIdTime], distAmongDrones
	%%
	% yyaxis left;
	% boxplot(squeeze(experimentsStruct.(typeOfMotion).idTime(:,end,1,:)), confidenceThresholds, 'boxstyle','filled', 'medianstyle','target', 'sym','r.', 'positions',confidenceThresholds);
	% yyaxis right;
	% boxplot(squeeze(experimentsStruct.(typeOfMotion).idMOTA(:,end,1,:)), confidenceThresholds, 'boxstyle','filled', 'medianstyle','target', 'sym','r.', 'positions',confidenceThresholds);
	for N = 5:5:25
		Nind = N-1;
		for normTypeInd = 1 %:normalizationTypes
			for typeOfMotionCell = {'random', 'hovering'} % {'random', 'hovering', 'landed'}
				typeOfMotion = typeOfMotionCell{:};
				if strcmp(typeOfMotion, 'random'), tit = 'Random motion'; else, tit = typeOfMotion; end
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
		for typeOfMotionCell = {'random', 'hovering', 'landed'}
			typeOfMotion = typeOfMotionCell{:};
			if strcmp(typeOfMotion, 'random'), tit = 'Random motion'; else, tit = typeOfMotion; end
			figure;
			bar(squeeze(experimentsStruct.(typeOfMotion).distAmongDrones(1,Nind,normTypeInd,:,1)), squeeze(experimentsStruct.(typeOfMotion).distAmongDrones(1,Nind,normTypeInd,:,2)));
			title([tit ', N=' num2str(N)], 'FontSize',18); xlabel('Distance among drones (m)', 'FontSize',18); ylabel('Count (#)', 'FontSize',18); set(gca, 'FontSize',14);
			saveas(gcf, ['histDistAmongDrones_' typeOfMotion '_N' num2str(N)], 'epsc');
		end
	end
end