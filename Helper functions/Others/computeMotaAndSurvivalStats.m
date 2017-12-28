function [MOTAmean, MOTAstd, numSurvivorsMean, numSurvivorsStd] = computeMotaAndSurvivalStats(experimentsStruct, typeOfMotion)
	MOTAmean = 100*reshape(mean(experimentsStruct.(typeOfMotion).percentCorrectOverTime(:,:,1,:), 1, 'omitnan'), size(experimentsStruct.(typeOfMotion).idMOTA,2), []);
	MOTAstd = 100*reshape(std(experimentsStruct.(typeOfMotion).percentCorrectOverTime(:,:,1,:), 0,1, 'omitnan'), size(experimentsStruct.(typeOfMotion).idMOTA,2), []);
	numSurvivorsMean = reshape(mean(experimentsStruct.(typeOfMotion).numSurvivors, 1, 'omitnan'), size(experimentsStruct.(typeOfMotion).numSurvivors,2), []);
	numSurvivorsStd = reshape(std(experimentsStruct.(typeOfMotion).numSurvivors, 0,1, 'omitnan'), size(experimentsStruct.(typeOfMotion).numSurvivors,2), []);
end
