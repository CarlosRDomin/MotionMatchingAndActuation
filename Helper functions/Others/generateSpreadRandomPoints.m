function pos = generateSpreadRandomPoints(N, minD, roomDimensions, wallMargin)
	if nargin<4 || isempty(wallMargin), wallMargin = 1; end

	while 1  % Keep trying until we get a valid distribution
		pos = generateValidPoint(minD, roomDimensions, wallMargin);
		for nTries = 1:20*N	% Limit the number of tries so we don't get stuck in a bad combination
			potentialNewPoint = generateValidPoint(minD, roomDimensions, wallMargin);
			if min(pdist2(potentialNewPoint, pos)) > minD
				pos = [pos; potentialNewPoint];
				if size(pos,1) >= N, return; end
			end
		end
	end
end

function point = generateValidPoint(minD, roomDimensions, wallMargin)
	if nargin<3 || isempty(wallMargin), wallMargin = 1; end

	point = wallMargin*minD + rand(1,2).*(roomDimensions(1:2)-2*wallMargin*minD);
end
