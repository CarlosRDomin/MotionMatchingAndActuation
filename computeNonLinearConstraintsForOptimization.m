%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the non-linear risk constraints (ie, drone-drone collisions) for a command given an assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [constrSmallerThan0, constrEqualTo0 , gradConstrSmallerThan0, gradConstrEqualTo0] = computeNonLinearConstraintsForOptimization(command, assignments, posUAV, threshRisk)
	newPosUAV = posUAV;	% Mx3 matrix
	newPosUAV(assignments(:,2),:) = posUAV(assignments(:,2),:) + reshape(command, 3,[])';
	N = size(posUAV,1);
	Ndist = N*(N-1)/2;

	cnt = 0;
	newPosSet1 = zeros(Ndist, 3);
	newPosSet2 = zeros(Ndist, 3);
	gradConstrSmallerThan0 = zeros(3*N, Ndist);
	gradAuxMatrix = zeros(3, 3*N, Ndist);
	for i = 1:N
		j = N-i;
		% The next j entries will be the dist between drone i and drones i+1:N
		newPosSet1(cnt + (1:j),:) = repmat(newPosUAV(i,:), j,1);
		newPosSet2(cnt + (1:j),:) = newPosUAV(i+(1:j),:);
		% And the aux matrix for the gradient should have -2s in the diagonal corresponding to command_i, and 2s in the diagonals corresponding to each command_j
		gradAuxMatrix(:, 3*(i-1) + (1:3), cnt + (1:j)) = repmat(-2*eye(3), 1,1,j);
		for k = 1:j
			gradAuxMatrix(:, 3*(i+k-1) + (1:3), cnt + k) = 2*eye(3);
		end
		% Finally, increase the counter so we don't overwrite data
		cnt = cnt+j;
	end
	newPdistPerAxis = newPosSet1 - newPosSet2;
	constrSmallerThan0 = threshRisk.^2 - sum(newPdistPerAxis.^2, 2);
	for i = 1:Ndist
		gradConstrSmallerThan0(:,i) = newPdistPerAxis(i,:)*gradAuxMatrix(:,:,i);
	end
	
% 	[~, distAmongDrones, distToWalls, newPosUAV] = estimateRiskOfCommand(command, assignments, posUAV, roomDimensions, false);
% 	constrSmallerThan0 = threshRisk - [distAmongDrones(:)]; %; distToWalls(:)];
	constrEqualTo0 = [];
	gradConstrEqualTo0 = [];
end
