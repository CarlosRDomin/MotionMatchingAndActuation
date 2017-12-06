function pos = generateSpreadRandomPoints(N, minD, roomDimensions)
	pos = generateValidPoint(minD, roomDimensions);
	k = 1;
	while (k < N)
		rv = generateValidPoint(minD, roomDimensions);
		if min(pdist2(rv, pos)) > minD
			pos = [pos; rv];
			k = k+1;
		end
	end
end

function point = generateValidPoint(minD, roomDimensions, wallMargin)
	if nargin<3 || isempty(wallMargin), wallMargin = 0.8; end
	point = wallMargin*minD + rand(1,2).*(roomDimensions(1:2)-2*wallMargin*minD);
end
