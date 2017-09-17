%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs our matching framework over a set of cam+IMU sensed accelerations by repeatedly calling computeBayesianIteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = runMatchingFrameworkOnGivenData(data, iUAVs, iCams, runningCorrWinSizes, dims, cropToMinT, maxExperimentTime)
	if nargin<2 || isempty(iUAVs)
		iUAVs = 1:length(data);
	end
	if nargin<3 || isempty(iCams)
		iCams = 1:length(data);
	end
	if nargin<4 || isempty(runningCorrWinSizes)
		runningCorrWinSizes = [25, 50, 75, 100, 250]; % [15:15:300 10000];
	end
	if nargin<5 || isempty(dims)
		dims = 1:3;
	end
	if nargin<6 || isempty(cropToMinT)
		cropToMinT = false;
	end
	if nargin<7 || isempty(maxExperimentTime)
		maxExperimentTime = NaN;
	end
	
	N = length(iUAVs);	% Number of drones transmitting IMU data
	M = length(iCams);	% Number of objects seen by a spotter
	minT = NaN; maxT = NaN;
	for iUAV = 1:length(data)
		tUAV = data(iUAV).a_UAV.X.t';
		tCam = data(iUAV).a_cam.X.t';
		[minT, indMinT] = min([minT, tUAV(end), tCam(end)]);
		[maxT, indMaxT] = max([maxT, tUAV(end), tCam(end)]);
		if indMinT==2, tMinT=tUAV; elseif indMinT==3, tMinT=tCam; end
		if indMaxT==2, tMaxT=tUAV; elseif indMaxT==3, tMaxT=tCam; end
	end
	if cropToMinT, t = tMinT; else, t = tMaxT; end	% If requested, crop all signals to the shortest signal. Otherwise, zero-pad all signals to the longest signal
	if ~isnan(maxExperimentTime), t(t>maxExperimentTime)=[]; end	% If requested, limit the maximum experiment time
	
	outFields = {'runningWinScore','runningLikelihood','runningPrior','assignedMatch','N','M','iCams','dims','runningCorrWinSizes','cropToMinT','t','yCam','yUAV'};
	out = cell2struct(cell(1, length(outFields)), outFields, 2);
	%%% t = zeros(M, lenT, length(dims));
	yCam = NaN(M, length(t), length(dims));
	yUAV = NaN(N, length(t), length(dims));
	runningWinScore = NaN(N, M, length(t), length(dims), length(runningCorrWinSizes));
	runningLikelihood = NaN(N, M, length(t), length(runningCorrWinSizes));
	runningPrior = cat(3, ones(N, M, 2, length(runningCorrWinSizes))./(N+M-1), NaN(N, M, length(t)-1, length(runningCorrWinSizes)));
	assignedMatch = NaN(N, length(t), length(runningCorrWinSizes));
	
	for iD = dims
		d = dims(iD);
		strAx = char('X'+d-1);	% Letter representation of the dimension ('X', 'Y', 'Z')
		for iC = 1:M
			iCam = iCams(iC);
			yCam(iC,:,iD) = interp1(data(iCam).a_cam.(strAx).t, data(iCam).a_cam.(strAx).measured, t);
		end
		for iU = 1:N
			iUAV = iUAVs(iU);
			yUAV(iU,:,iD) = interp1(data(iUAV).a_UAV.(strAx).t, data(iUAV).a_UAV.(strAx).measured, t);
		end
	end

	dispImproved('', 'init'); dispImproved('Computing stats... ', 'keepthis');
	for currT = 2:length(t)
		[runningWinScore, runningLikelihood, runningPrior, assignedMatch] = computeBayesianIteration(runningWinScore, runningLikelihood, runningPrior, assignedMatch, yCam, yUAV, currT, dims, runningCorrWinSizes, N, M);
		dispImproved(sprintf('t=%6.2fs; %6.2f%% (%d out of %d)\n', t(currT), 100*currT/length(t), currT, length(t)));
	end
	for f = outFields	% Populate output struct with results
		out.(f{:}) = eval(f{:});
	end
end
