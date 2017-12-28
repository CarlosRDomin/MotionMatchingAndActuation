%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the safest direction for all drones to move at once (cross product of the 2 shortest drone pairs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rhoBounds, dirTheta, dirPhi] = estimateLowestRiskyDirOfMotion(posUAV, roomDimensions, threshRisk, deltaP)
	% Compute the distance among pairs of drones and from each drone to each wall
	distAmongDrones = pdist(posUAV);
	distToWalls = [posUAV - zeros(size(posUAV)), repmat(roomDimensions, size(posUAV,1),1) - posUAV];
	
	% Then combine all distances and sort them (to find minimum distances)
	d = sort([distAmongDrones(:); distToWalls(:)]); %d(d<=0) = []; % Make sure distances are valid
	squareDistAmongDrones = squareform(distAmongDrones); % squareform() converts vector to matrix -> We can easily know indices of shortest distance pair
	triUppDistAmongDrones = tril(NaN(size(squareDistAmongDrones))) + squareDistAmongDrones; % tril() keeps the lower triangular part (so indices aren't found twice [i,j], [j,i])
	
	shortestVectors = NaN(2,3+1); % Each row contains a 3d vector with a specific direction + a '0' if it corresponds to a drone-drone case, '1' if a wall
	shortestVectorRowInd = 1;
	iD = 1;
	while ~testShortestVectorsValid(shortestVectors) % Keep finding shortest vectors until the cross product of them two is perpendicular to them (avoid 2 parallel shortest vectors)
		shouldBreak = false;
		% Test to see if any pair-wise drone distance matches d(iD) [otherwise d(iD) is the distance from one drone to one of the walls, will test that next]
		[indAmongDrones1, indAmongDrones2] = find(triUppDistAmongDrones==d(iD));
		for iFound = 1:length(indAmongDrones1) % In theory, there could be multiple pair-wise distances all equal to d(iD). Iterate through them
			shortestVectors(shortestVectorRowInd,:) = [posUAV(indAmongDrones1(iFound),:)-posUAV(indAmongDrones2(iFound),:), 0];
			if testShortestVectorsValid(shortestVectors), shouldBreak=true; break; end
			shortestVectorRowInd = 2;
		end
		if shouldBreak, break; end % The inner for-loop signals that shortestVectors is valid and we can exit this while-loop
		
		% Test to see if any of the drones is at distance d(iD) to any of the walls
		[~, indWall] = find(distToWalls==d(iD));
		for iFound = 1:length(indWall) % In theory, there could be multiple drone-to-wall distances all equal to d(iD). Iterate through them
			whichAxis = (1:3 == mod(indWall(iFound)-1,3)+1); % indWall is a 1:6 number that determines the shortestVector axis directly (eg, 1 and 4 are x-axis [YZ plane])
			whichAxisWithSign = ((-1)^(indWall(iFound)>3)) * whichAxis; % whichAxisWithSign points towards inside of the room -> Eg: wall 1 is [1,0,0] but wall 4 is [-1,0,0]
			shortestVectors(shortestVectorRowInd,:) = [max(d(iD), 0.1).*whichAxisWithSign, 1]; % Scale whichAxisWithSign by the distance in case we have to linearly combine 2 vectors
			if testShortestVectorsValid(shortestVectors), shouldBreak=true; break; end
			shortestVectorRowInd = 2;
		end
		if shouldBreak, break; end % The inner for-loop signals that shortestVectors is valid and we can exit this while-loop
		
		% Keep iterating the d vector
		iD = iD + 1;
		if iD > length(d), error('This should never happen: unable to estimate lowest risky direction of motion :('); end
	end
	
	%maxRho = d(1); % This is very safe, should be the shortest distance in the direction of motion
	dirXYZ = computeDirXYZ(shortestVectors);
	[dirXYZ, rhoBounds] = findBestWayInDirection(dirXYZ, posUAV, triUppDistAmongDrones, distToWalls, roomDimensions, threshRisk, [], deltaP);
	[dirTheta,dirPhi] = cart2sph(dirXYZ(1), dirXYZ(2), dirXYZ(3));
	dirPhi = pi/2 - dirPhi;
end

function dirXYZ = computeDirXYZ(shortestVectors)
	vectors = shortestVectors(:,1:3);
	vectorTypes = shortestVectors(:,4);
	
	if all(vectorTypes == 0)		% 2 drone-drone distances: resulting direction is perpendicular to both
		dirXYZ = cross(vectors(1,:), vectors(2,:));
	elseif all(vectorTypes == 1)	% 2 walls: resulting direction is linear combination of both vectors (away from the wall)
		norm1 = norm(vectors(1,:));
		norm2 = norm(vectors(2,:));
		dirXYZ = norm2/norm1*vectors(1,:) + norm1/norm2*vectors(2,:);
	else % 1 drone-drone, 1 wall
		dirXYZ = projectVectPlane(vectors(vectorTypes==1,:), vectors(vectorTypes==0,:));
	end
end

function isValid = testShortestVectorsValid(shortestVectors)
	isValid = (all(~isnan(shortestVectors(:,4))) && norm(computeDirXYZ(shortestVectors)) > 1e-5);
end

function [dirXYZ, rhoBounds] = findBestWayInDirection(dirXYZ, posUAV, triUppDistAmongDrones, distToWalls, roomDimensions, threshRisk, riskMode, deltaP)
	if nargin<7 || isempty(riskMode), riskMode = 1; end % riskMode: 0 is very conservative (when in doubt, hover), 1 is less conservative (drone-drone constraints already within threshRisk are ignored), 2 is way less conservative (maxRho is the rho that would hit a drone)
	if nargin<8, deltaP = []; end
	if riskMode==0, distDronesAlreadyInsideThreshRisk = 0; elseif riskMode==1, distDronesAlreadyInsideThreshRisk = NaN; end
	
	% Find the direction (not the "way") as the cross product of the 2 shortest vectors
	dirXYZ = dirXYZ./norm(dirXYZ); % Make unitary
	rhoBounds = NaN(1,2); % maxRho(1) is in the positive +dirXYZ way; maxRho(2) is in -dirXYZ
	
	% Find (if any) restricted ways: if there's a drone closer to a wall than threshRisk, set the corresponding maxRho to 0 to avoid it getting even closer to the wall
	angleOfDirXYZwrtWalls = sign(dot(repmat(dirXYZ, 6,1),[eye(3);-eye(3)], 2)'); % Use the sign of the dot product dirXYZ•wallAwayVect to determine whether +dirXYZ brings you closer or away from the wall
	for riskyWallsInd = find(any(distToWalls<threshRisk, 1)) % For every wall with drones closer than threshRisk
		if angleOfDirXYZwrtWalls(riskyWallsInd) ~= 0 % Corner case, if dirXYZ is perpendicular to the wall, then no need to take action (we're not gonna get closer)
			rhoBounds(1+(angleOfDirXYZwrtWalls(riskyWallsInd)==1)) = 0; % If dirXYZ takes you away from the wall (positive dot product), maxRho in the -dirXYZ direction (maxRho(2)) should be 0; and viceversa: if dirXYZ takes you into the wall (negative dot product), do not allow to go in +dirXYZ direction (maxRho(1)=0)
		end
	end
	
	% Now, we need to figure out which way to go to maximize maxRho.
	lambdaDrones = []; % Initialize in case riskMode is 2 (= ignore lambdaDrones)
	planeVectors = [0 1 0, 0 0 1; 1 0 0, 0 0 1; 1 0 0, 0 1 0]; % x-axis plane = YZ plane -> vectors are y-axis and z-axis; y-axis plane = XZ; z-axis plane = XY
	for i = 1:size(posUAV,1)
		% The shortest distance is the one where one drone doesn't move at all, and the other one moves the maximum allowed.
		% We want that distance to be at least threshRisk, so intersect the line representing drone i moving in dirXYZ with the sphere at all other drones i+1:end and Radius=threshRisk
		lineDroneI = [posUAV(i,:) dirXYZ];
		if riskMode < 2
			[~,lambdaDrones] = intersectLineSphere(lineDroneI, [posUAV(i+1:end,:) repmat(threshRisk, size(posUAV,1)-i,1)]);
			lambdaDrones(reshape(repmat(triUppDistAmongDrones(i,i+1:end), 2,1), [],1) < threshRisk) = distDronesAlreadyInsideThreshRisk; % Decide how to react if both drones are already within threshRisk
		end
		[~,lambdaWalls] = intersectLinePlane(lineDroneI, [repmat(threshRisk, 3,3) planeVectors; repmat(roomDimensions-threshRisk, 3,1) planeVectors]);
		lambdaWalls(distToWalls(i,:)<threshRisk) = 0; % Ignore lambdaWalls corresponding to drones outside the "safe zone" (already took care of it outside this for-loop)
		lambda = [lambdaDrones(:); lambdaWalls(:)]; % Combine results found for both drone-drone and drone-wall collisions
		
		% Now, find the minimum between the current minimum distance in each way, and the new constraints obtained for this drone
		minPosLambda = min(lambda(sign(lambda)==1));
		if ~isempty(minPosLambda), rhoBounds(1) = min(rhoBounds(1), minPosLambda); end
		minNegLambda = min(abs(lambda(sign(lambda)==-1)));
		if ~isempty(minNegLambda), rhoBounds(2) = min(rhoBounds(2), minNegLambda); end
	end
	
	% Finally, choose which direction is less risky (higher maxRho)
	rhoBounds(isnan(rhoBounds))=0; % Testing for > or < against a NaN always returns false. Fix that by replacing NaN by 0.
	if rhoBounds(2) > rhoBounds(1)
		rhoBounds = rhoBounds(2:-1:1); % Let's always have rhoBounds(1) be greater than rhoBounds so dirXYZ indicates dir of max motion
		dirXYZ = -dirXYZ;
	end
	% rhoBounds = sign(rhoBounds).*max(abs(rhoBounds), threshRisk/4);
	if ~isempty(deltaP), rhoBounds = min(rhoBounds, deltaP); end % Make sure it's a feasible command
end
