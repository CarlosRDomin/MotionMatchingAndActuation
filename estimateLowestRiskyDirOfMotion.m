%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the safest direction for all drones to move at once (cross product of the 2 shortest drone pairs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [maxRho, dirTheta, dirPhi] = estimateLowestRiskyDirOfMotion(posUAV)
	d = pdist(posUAV);
	[~,iShortest] = sort(d);
	D = squareform(d);	% squareform() converts vector to matrix -> We can easily know indices of shortest distance pair
	[i1,j1] = find(D==d(iShortest(1)), 1);
	for i = 2:length(d)
		[i2,j2] = find(D==d(iShortest(i)), 1);
		dirXYZ = cross(posUAV(i1,:)-posUAV(j1,:), posUAV(i2,:)-posUAV(j2,:));
		if norm(dirXYZ)>1e-3, break; end;
	end
	
	[dirTheta,dirPhi] = cart2sph(dirXYZ(1),dirXYZ(2),dirXYZ(3));
	dirPhi = pi/2 - dirPhi;
	maxRho = d(iShortest(1));
end
