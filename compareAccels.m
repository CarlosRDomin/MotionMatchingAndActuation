%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates different plots from a set of real experiments (simply calls plotExperimentResults)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
cdToThisScriptsDirectory;

runningCorrStruct = plotExperimentResults([], struct('rawAccel',false, 'scatterCorr',false, 'runningLikelihoodVsWinSize',false, 'runningLikelihoodFull',true));