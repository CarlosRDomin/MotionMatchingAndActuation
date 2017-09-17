%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicts where drones would end up after a given command and assignment and returns the new shortest distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function risk = estimateRiskOfCommand(command, assignments, posUAV)
	newPosUAV = posUAV;	% Mx3 matrix
	[dx,dy,dz] = sph2cart(command(3*assignments(:,2)-1), pi/2-command(3*assignments(:,2)), command(3*assignments(:,2)-2));
	newPosUAV(assignments(:,2),:) = posUAV(assignments(:,2),:) + [dx,dy,dz];
	d = pdist(newPosUAV);
	risk = min(d);
end
