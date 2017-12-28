function command = deltaPtoCommand(deltaP)
	% Compute which command would lead to each drone moving by deltaP=[dx dy dz] -> N x 3
	[commandTheta, commandElev, commandRho] = cart2sph(deltaP(:,1), deltaP(:,2), deltaP(:,3));
	commandPhi = pi/2-commandElev;
	command = reshape([commandRho'; commandTheta'; commandPhi'], [],1);
end
