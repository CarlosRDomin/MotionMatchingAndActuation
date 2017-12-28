function deltaP = commandToDeltaP(commandRho, commandTheta, commandPhi)
	% Compute how much we want to move in each cartesian direction
	commandRho = reshape(commandRho, [],1);
	commandTheta = reshape(commandTheta, [],1);
	commandElev = reshape(pi/2 - commandPhi, [],1);
	[dx,dy,dz] = sph2cart(commandTheta, commandElev, commandRho); % Convert spherical to cartesian coords
	deltaP = [dx dy dz]; % numel(commandRho) x 3
end
