function out = computeRunningCorr(data, iCams, runningCorrWinSizes, dims, cropToMinT)
	if nargin<2 || isempty(iCams)
		iCams = 1:length(data);
	end
	if nargin<3 || isempty(runningCorrWinSizes)
		runningCorrWinSizes = [25, 50, 75, 100, 250]; % [15:15:300 10000];
	end
	if nargin<4 || isempty(dims)
		dims = 1:3;
	end
	if nargin<5 || isempty(cropToMinT)
		cropToMinT = false;
	end
	
	N = length(data);	% Number of drones transmitting IMU data
	M = length(iCams);	% Number of objects seen by a spotter
	minT = NaN; maxT = NaN;
	for iUAV = 1:length(data)
		tUAV = data(iUAV).a_UAV.X.t';
		tCam = data(iUAV).a_cam.X.t';
		[minT, indMinT] = min([minT, tUAV(end), tCam(end)]);
		[maxT, indMaxT] = max([maxT, tUAV(end), tCam(end)]);
		if indMinT==2, tMinT=tUAV; else, tMinT=tCam; end
		if indMaxT==2, tMaxT=tUAV; else, tMaxT=tCam; end
	end
	if cropToMinT, t = tMinT; else, t = tMaxT; end	% If requested, crop all signals to the shortest signal. Otherwise, zero-pad all signals to the longest signal
	t(t>30)=[];
	
	outFields = {'runningCorr','runningLikelihood','runningPrior','assignedMatch','N','M','iCams','dims','runningCorrWinSizes','cropToMinT','t','yCam','yUAV'};
	out = cell2struct(cell(1, length(outFields)), outFields, 2);
	%%% t = zeros(M, lenT, length(dims));
	yCam = NaN(M, length(t), length(dims));
	yUAV = NaN(N, length(t), length(dims));
	runningCorr = NaN(N, M, length(t), length(dims), length(runningCorrWinSizes));
	assignedMatch = NaN(M, length(t), length(dims), length(runningCorrWinSizes));
	muSgivenTheta = eye(N); sigmaSgivenTheta = repmat(1.5*eye(N), 1,1,N);
	runningPrior = cat(3, ones(N, M, 2, length(dims), length(runningCorrWinSizes))./N, NaN(N, M, length(t)-1, length(dims), length(runningCorrWinSizes)));
	runningLikelihood = NaN(N, M, length(t), length(dims), length(runningCorrWinSizes));
	
	for iD = dims
		d = dims(iD);
		strAx = char('X'+d-1);	% Letter representation of the dimension ('X', 'Y', 'Z')
		for iC = 1:M
			iCam = iCams(iC);
			yCam(iC,:,iD) = interp1(data(iCam).a_cam.(strAx).t, data(iCam).a_cam.(strAx).measured, t);
		end
		for iUAV = 1:N
			yUAV(iUAV,:,iD) = interp1(data(iUAV).a_UAV.(strAx).t, data(iUAV).a_UAV.(strAx).measured, t);
		end
	end

	dispImproved('', 'init'); dispImproved('Computing stats... ', 'keepthis');
	for currT = 2:length(t)
		[runningCorr, runningLikelihood, runningPrior, assignedMatch] = computeBayesianIteration(runningCorr, runningLikelihood, runningPrior, assignedMatch, yCam, yUAV, currT, dims, runningCorrWinSizes, N, M, muSgivenTheta, sigmaSgivenTheta);
		dispImproved(sprintf('t=%6.2fs; %6.2f%% (%d out of %d)\n', t(currT), 100*currT/length(t), currT, length(t)));
	end
	for f = outFields	% Populate output struct with results
		out.(f{:}) = eval(f{:});
	end
end
