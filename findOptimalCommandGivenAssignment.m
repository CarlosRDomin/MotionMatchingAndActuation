%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the best command to give each drone so that the uniqueness in their motion is maximized, subject to risk constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [command,newP,exitFlag,output] = findOptimalCommandGivenAssignment(assignments, posUAV, roomDimensions, threshRisk, runningPrior, spotterCam, deltaP, deltaT, sigmaNoiseAccel, startingCommand)
	if ndims(posUAV) == 3, posUAV = reshape(posUAV, [],size(posUAV,3)); end % Should already be Mxlength(dims) but sometimes it's Mx1xlength(dims) (because we have an array of Mxlength(t)xlength(dims))
	N = size(runningPrior, 1);

	% Linear constraints: drone-wall collisions
	% Command vector x has length 3*N: c1_x, c1_y, c1_z, c2_x, c2_y, c2_z, c3_x...
	% There's 2 walls per x/y/z dimension -> 2*3*N linear constraints. They can be divided in 2 sets:
	%  - "Min" walls (W_xMin, W_yMin, W_zMin) -> Must satisfy: W_dMin + D_th <= Pi_d + Ci_d  ->  -Ci_d <= Pi_d - W_dMin - D_th
	%  - "Max" walls (W_xMax, W_yMax, W_zMax) -> Must satisfy: Pi_d + Ci_d <= W_dMax - D_th  ->   Ci_d <= W_dMax - D_th - Pi_d
	posUAVcolumn = reshape(posUAV', [],1); % posUAV unraveled: p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, p3_x...
	WmaxMinusDth = repmat(reshape(roomDimensions-threshRisk, [],1), N,1); % W_dMax-D_th unraveled: W_xMax-D_th, W_yMax-D_th, W_zMax-D_th, W_xMax-D_th...
	
	A = [-eye(3*N); eye(3*N)]; % W_min constraints are 0 0 ...0 -1 0... 0 0, one per command element; W_max are 0 0 ...0 +1 0... 0 0
	b = [posUAVcolumn-threshRisk; WmaxMinusDth-posUAVcolumn];

	% Use Matlab's optimization solver
	[command,newP,exitFlag,output] = fmincon(@(x) (-estimateImprovementOfCommand(deltaPtoCommand(reshape(x, 3,[])'),assignments,runningPrior,[],sigmaNoiseAccel,[],spotterCam.fps,deltaT)), ...
		zeros(3*N,1),A,b,[],[],-deltaP.*ones(3*N,1),deltaP.*ones(3*N,1), ...
		@(x) computeNonLinearConstraintsForOptimization(x, assignments, posUAV, threshRisk), optimoptions('fmincon', 'FunctionTolerance',1e-3, 'StepTolerance',1e-3, 'OptimalityTolerance',1e-8, 'SpecifyConstraintGradient',true)); %, 'StepTolerance',1e-3, 'Algorithm','sqp', 'FunctionTolerance',1e-3, 'UseParallel',true, 'Display','iter-detailed'));
	
	% Old code:
	% [command,newP,exitFlag,output] = fmincon(@(x) (-estimateImprovementOfCommand(x,assignments,runningPrior,[],sigmaNoiseAccel,[],spotterCam.fps,deltaT)), ...
	% 	startingCommand,[],[],[],[],zeros(3*N,1),repmat([deltaP,2*pi,pi]',N,1), ...
	% 	@(x) getNonLinearConstraintsForOptimization(x, assignments, posUAVcam, roomDimensions, threshRisk), optimoptions('fmincon', 'FunctionTolerance',1e-3)); %, 'Algorithm','sqp', 'FunctionTolerance',1e-3, 'UseParallel',true, 'Display','iter-detailed'));
end
