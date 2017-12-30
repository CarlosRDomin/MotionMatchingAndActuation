%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicts where drones would end up after a given command and assignment and returns the new shortest distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [risk, distAmongDrones, distToWalls, newPosUAV] = estimateRiskOfCommand(command, assignments, posUAV, roomDimensions, isCommandSpherical)
	if nargin<5 || isempty(isCommandSpherical), isCommandSpherical = true; end

	commandMat = reshape(command, 3,[])'; % Convert to Nx3 matrix form so we can easily access all rhos, all thetas, etc.
	if ndims(posUAV) == 3, posUAV = reshape(posUAV, [],size(posUAV,3)); end % Should already be Mxlength(dims) but sometimes it's Mx1xlength(dims) (because we have an array of Mxlength(t)xlength(dims))
	newPosUAV = posUAV;	% Mx3 matrix
	if isCommandSpherical
		[dx,dy,dz] = sph2cart(commandMat(:,2), pi/2-commandMat(:,3), commandMat(:,1));
	else
		dx = commandMat(:,1); dy = commandMat(:,2); dz = commandMat(:,3);
	end
	newPosUAV(assignments(:,2),:) = posUAV(assignments(:,2),:) + [dx,dy,dz];
	distAmongDrones = pdist(newPosUAV);
	distToWalls = [newPosUAV - zeros(size(newPosUAV)), repmat(roomDimensions, size(newPosUAV,1),1) - newPosUAV];
	risk = min([distAmongDrones(:); distToWalls(:)]);
end
