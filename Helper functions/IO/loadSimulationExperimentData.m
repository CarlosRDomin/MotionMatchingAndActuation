function [expData, experimentInd] = loadSimulationExperimentData(experimentStruct)
	if isfield(experimentStruct, 'experimentMatData') % If the data is already avaialble, simply use the data directly
		expData = experimentStruct.experimentMatData;
	else % Otherwise, load the *.mat
		expData = load(experimentStruct.experimentMatFileName);
	end
	experimentInd = experimentStruct.experimentInd; % Keep track of which experiment within all repsPerExperiment rows we want to plot
end
