%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots all figures related to the matching framework evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; %clear all;
cdToThisScriptsDirectory;

simOutputFolder = [fileparts(mfilename('fullpath')) '/data/Simulations/'];
typeOfExperiment = 'MatchingFramework';
filesInFolder = dir([simOutputFolder typeOfExperiment '_*.mat']);

experimentsStruct = struct();
for fileNameCell = {filesInFolder.name}
	fileName = fileNameCell{:};
	fileInfo = strsplit(fileName, '_');
	if strcmp(fileInfo{2}, 'processed.mat'), continue; end % Skip processed results
	N = str2double(fileInfo{3}(1:end-1));
	if strcmp(fileInfo{1}, typeOfExperiment) && N>1 && N<21 % Only process files of the right type of experiment
		if ~isfield(experimentsStruct, fileInfo{2}) % Then, arrange the *.mat's according to their motion type (fileInfo{2}, ie. random, hovering, landed...)
			experimentsStruct.(fileInfo{2}).fileNames = {fileName};	     % Initialize the field if it didn't exist
			experimentsStruct.(fileInfo{2}).N = N;
		else
			experimentsStruct.(fileInfo{2}).fileNames{end+1} = fileName; % Otherwise just append the fileName to the end of the list
			experimentsStruct.(fileInfo{2}).N(end+1) = N;
		end
	end
end

%%
confidenceThresholds = 0:0.01:1;
repsPerExperiment = 20;
normalizationTypes = 3;
distHistNbins = 50;

for typeOfMotionCell = fieldnames(experimentsStruct)'
	typeOfMotion = typeOfMotionCell{:};
	
	% Sort fileNames by N, and save N as a 
	[experimentsStruct.(typeOfMotion).N, sortedInds] = sort(experimentsStruct.(typeOfMotion).N);
	experimentsStruct.(typeOfMotion).N = unique(experimentsStruct.(typeOfMotion).N); % There might be repeated values of N (different normalizations), so keep only unique values
	experimentsStruct.(typeOfMotion).fileNames = experimentsStruct.(typeOfMotion).fileNames(sortedInds); % Sort fileNames by N

	% Initialize remaining fields in the structure
	experimentsStruct.(typeOfMotion).idAccuracy = NaN(repsPerExperiment, length(experimentsStruct.(typeOfMotion).N), normalizationTypes, length(confidenceThresholds));
	experimentsStruct.(typeOfMotion).idTime = NaN(size(experimentsStruct.(typeOfMotion).idAccuracy));
	experimentsStruct.(typeOfMotion).distAmongDrones = zeros(repsPerExperiment, length(experimentsStruct.(typeOfMotion).N), normalizationTypes, distHistNbins, 2);

	for fileNameInd = 1:numel(experimentsStruct.(typeOfMotion).fileNames)
		fileName = experimentsStruct.(typeOfMotion).fileNames{fileNameInd};
		expData = load(fileName); % Load the experiment mat
		N = expData.paramStruct.N; [~,Nind] = find(experimentsStruct.(typeOfMotion).N==N);
		fileInfo = strsplit(fileName, '_'); normType = str2double(fileInfo{4}(end)); normTypeInd = normType+1; % Forgot to save normType into the experiment parameters, so extract it from file name ;)
		
		% Traverse every experiment repetition for these variables
		for experimentInd = 1:length(expData.variableStruct)
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
				correctOverTime = sum(confidentAssignments==gtAssignment, 1, 'omitNaN');
				percentCorrectOverTime = correctOverTime ./ max(1, sum(~isnan(confidentAssignments), 1)); % Make sure we don't divide 0/0. Make divisor >=1 so output would be 0/1=0 in that case
				experimentsStruct.(typeOfMotion).idAccuracy(experimentInd,Nind,normTypeInd,confidenceThreshInd) = sum(correctOverTime(:)) / max(1, sum(~isnan(confidentAssignments(:))));
				indFirstAllIDd = find(correctOverTime == N, 1); % It could be possible that (eg, for high confidenceThresh values) the system never fully IDs all drones. idTime is initialized with NaNs so just don't write to it in that case
				if ~isempty(indFirstAllIDd), experimentsStruct.(typeOfMotion).idTime(experimentInd,Nind,normTypeInd,confidenceThreshInd) = expData.variableStruct(experimentInd).t(indFirstAllIDd); end
				dispImproved(sprintf('\nProcessing experiment results for motion "%s":\tN=%2d, norm=%d (%2d/%2d)\texperiment=%2d/%2d\tconfidenceThresh=%3d%%...', typeOfMotion, N, normType, fileNameInd, numel(experimentsStruct.(typeOfMotion).fileNames), experimentInd, length(expData.variableStruct), round(100*confidenceThresh)));
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
processedResultsFileName = [simOutputFolder typeOfExperiment '_processed.mat'];
save(processedResultsFileName, 'experimentsStruct', '-v7.3');
dispImproved(sprintf('Saved processed experiment results as "%s"!\n', processedResultsFileName), 'keepthis');

%%
paramFields = {'N','M','dims','runningCorrWinSizes','iCams'};
variableFields = {'runningWinScore','runningLikelihood','runningPosterior','assignedMatch','t','accelCam','accelUAV'};
outFields = cat(2, paramFields, variableFields);
runningCorrStruct = cell2struct(cell(1, length(outFields)), outFields, 2);
runningCorrStruct.N = expData.paramStruct.N;
runningCorrStruct.M = runningCorrStruct.N;
runningCorrStruct.iCams = expData.variableStruct(experimentInd).groundTruthAssignment;
runningCorrStruct.runningCorrWinSizes = expData.paramStruct.frameworkWinSize;
runningCorrStruct.dims = expData.paramStruct.dims;
for f = variableFields	% Populate output struct with results
	runningCorrStruct.(f{:}) = eval(['expData.variableStruct(experimentInd).' f{:}]);
end
plotExperimentResults(runningCorrStruct, struct('rawAccel',true, 'scatterCorr',false, 'runningLikelihoodVsWinSize',false, 'runningLikelihoodFull',true));

%%
generateDronesInRoomVideo([], expData.variableStruct(experimentInd).posUAVgt, expData.paramStruct.roomDimensions, expData.paramStruct.spotterCam, false, expData.paramStruct.frameworkWinSize);

%%
yyaxis left;
boxplot(squeeze(idTime(:,end,:)), confidenceThresholds, 'boxstyle','filled', 'medianstyle','target', 'sym','r.', 'positions',confidenceThresholds);
yyaxis right;
boxplot(squeeze(idAccuracy(:,end,:)), confidenceThresholds, 'boxstyle','filled', 'medianstyle','target', 'sym','r.', 'positions',confidenceThresholds);
