%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates different plots from a set of real experiments (simply calls plotExperimentResults)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
cdToWorkingDirectory;

runningCorrStruct = plotExperimentResults([], struct('rawAccel',false, 'scatterCorr',false, 'runningLikelihoodVsWinSize',true, 'runningLikelihoodFull',true));