%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicts where drones would end up after a given command and assignment and returns the new shortest distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [risk, distAmongDrones, distToWalls] = estimateRiskOfCommand(command, assignments, posUAV, roomDimensions)
	if ndims(posUAV) == 3, posUAV = reshape(posUAV, [],size(posUAV,3)); end % Should already be Mxlength(dims) but sometimes it's Mx1xlength(dims) (because we have an array of Mxlength(t)xlength(dims))
	newPosUAV = posUAV;	% Mx3 matrix
	[dx,dy,dz] = sph2cart(command(3*assignments(:,2)-1), pi/2-command(3*assignments(:,2)), command(3*assignments(:,2)-2));
	newPosUAV(assignments(:,2),:) = posUAV(assignments(:,2),:) + [dx,dy,dz];
	distAmongDrones = pdist(newPosUAV);
	distToWalls = [newPosUAV - zeros(size(newPosUAV)), repmat(roomDimensions, size(newPosUAV,1),1) - newPosUAV];
	risk = min(distAmongDrones); % min([distAmongDrones(:); distToWalls(:)]);
end
