function [constrSmallerThan0, constrEqualTo0] = getNonLinearConstraintsForOptimization(command, assignments, posUAVcam, roomDimensions, threshRisk)
	[~, distAmongDrones, distToWalls] = estimateRiskOfCommand(command, assignments, posUAVcam, roomDimensions);
	constrSmallerThan0 = threshRisk - [distAmongDrones(:); distToWalls(:)];
	constrEqualTo0 = 0;
end
