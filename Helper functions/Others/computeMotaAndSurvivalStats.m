function [MOTAmean, MOTAstd, numSurvivorsMean, numSurvivorsStd] = computeMotaAndSurvivalStats(experimentsStruct, typeOfMotion, intervalErrorBar)
	if nargin<3 || isempty(intervalErrorBar)
		intervalErrorBar = 2;
	end

	MOTAmean = 100*reshape(mean(experimentsStruct.(typeOfMotion).percentCorrectOverTime(:,:,1,:), 1, 'omitnan'), size(experimentsStruct.(typeOfMotion).idMOTA,2), []);
	MOTAstd = 100*reshape(std(experimentsStruct.(typeOfMotion).percentCorrectOverTime(:,:,1,:), 0,1, 'omitnan'), size(experimentsStruct.(typeOfMotion).idMOTA,2), []);
	numSurvivorsMean = 100*reshape(mean(experimentsStruct.(typeOfMotion).numSurvivors, 1, 'omitnan'), size(experimentsStruct.(typeOfMotion).numSurvivors,2), []);
	numSurvivorsStd = 100*reshape(std(experimentsStruct.(typeOfMotion).numSurvivors, 0,1, 'omitnan'), size(experimentsStruct.(typeOfMotion).numSurvivors,2), []);
	
	removeStdInds = 1:size(MOTAmean,2); removeStdInds(2:intervalErrorBar:end) = [];
	MOTAstd(:,removeStdInds) = 0; numSurvivorsStd(:,removeStdInds) = 0;
end
