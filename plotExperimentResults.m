%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates different plots (specified by the struct whatToPlot) from a set of real experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [runningCorrStruct] = plotExperimentResults(runningCorrStruct, whatToPlot, data)
	cdToThisScriptsDirectory;

	%% Load parameters / default values
	if nargin<1 || isempty(runningCorrStruct)
		try
			aux = load('runningCorrStruct.mat');
			runningCorrStruct = aux.runningCorrStruct;
			runningCorrStruct.frameworkWinSize([1:2 4:end]) = [];
			runningCorrStruct.runningWinScore(:,:,:,:,[1:2 4:end]) = [];
			runningCorrStruct.assignedMatch(:,:,:,[1:2 4:end]) = [];
			runningCorrStruct.runningPosterior(:,:,:,:,[1:2 4:end]) = [];
			runningCorrStruct.runningLikelihood(:,:,:,:,[1:2 4:end]) = [];
			runningCorrStruct.dims = 3;
			error;	% Force to go into loading experiment data from scratch
		catch
			if nargin<3 || isempty(data)	% Only load data if it hasn't been loaded yet (saves time)
				%data = loadRealExperimentData(struct('datetime',{'2017-02-19 17-14-15','2017-02-19 17-56-48','2017-02-19 17-59-23','2017-02-19 18-01-29','2017-02-19 18-22-23','2017-02-19 18-24-25'}, 'ch','75'));
				data = loadRealExperimentData(struct('datetime',{'2017-02-19 17-56-48','2017-02-19 17-59-23','2017-02-19 18-01-29','2017-02-19 18-22-23','2017-02-19 18-24-25'}, 'ch','75'));
			end
			runningCorrStruct = runMatchingFrameworkOnGivenData(data, [1 2 3 4 5], [1 2 3 4 5], [5:5:15 30:30:300]);	%(data, iCam=1:length(data), frameworkWinSize, dims=1:3, cropToMinT=true)
			save('runningCorrStruct.mat','runningCorrStruct');
		end
	elseif isfield(runningCorrStruct, 'experimentInd') && any(isfield(runningCorrStruct, {'experimentMatData', 'experimentMatFileName'})) % Allow to plot data directly from a simulation *.mat
		[expData, experimentInd] = loadSimulationExperimentData(runningCorrStruct);
		
		% Create a new runningCorrStruct and fill out the corresponding fields from the data loaded
		paramFields = {'N','M','dims','frameworkWinSize','iCams'};
		variableFields = {'runningWinScore','runningLikelihood','runningPosterior','assignedMatch','t','accelCam','accelUAV'};
		outFields = cat(2, paramFields, variableFields);
		runningCorrStruct = cell2struct(cell(1, length(outFields)), outFields, 2);
		runningCorrStruct.N = expData.paramStruct.N;
		runningCorrStruct.M = runningCorrStruct.N;
		runningCorrStruct.dims = expData.paramStruct.dims;
		runningCorrStruct.frameworkWinSize = expData.paramStruct.frameworkWinSize;
		runningCorrStruct.iCams = expData.variableStruct(experimentInd).groundTruthAssignment;
		for f = variableFields	% Populate output struct with results
			runningCorrStruct.(f{:}) = eval(['expData.variableStruct(experimentInd).' f{:}]);
		end
	end
	if nargin<2 || isempty(whatToPlot)
		whatToPlot = struct('rawAccel',true, 'scatterCorr',false, 'runningLikelihoodVsWinSize',true, 'runningLikelihoodFull',true);
	end

	figs = zeros(3, length(fieldnames(whatToPlot)));
	correctBgndColor = [0.9 1 0.9]; wrongBgndColor = [1 0.9 0.9];

	%% Plot rawAccel: plot in time-domain spotter's a_cam and UAV's a_UAV
	if whatToPlot.rawAccel
		figs(:, 1) = figure('Units','normalized', 'Position',[0 0 1 1]);
		ax = gobjects(runningCorrStruct.N, length(runningCorrStruct.dims));
		for iUAV = 1:runningCorrStruct.N		% For each IMU
			for iD = 1:length(runningCorrStruct.dims)
				d = runningCorrStruct.dims(iD);	% For each dimension (e.g. z)
				strAx = char('X'+d-1);			% Letter representation of the dimension ('X', 'Y', 'Z')
				ax(iUAV, iD) = subplot(runningCorrStruct.N,length(runningCorrStruct.dims), iD + (iUAV-1)*length(runningCorrStruct.dims)); hold on;
				if any(runningCorrStruct.iCams == iUAV)
					hCam = plot(runningCorrStruct.t, runningCorrStruct.accelCam(runningCorrStruct.iCams==iUAV,:,d), 'LineWidth',2);
				end
				hUAV = plot(runningCorrStruct.t, runningCorrStruct.accelUAV(iUAV,:,d)-mean(runningCorrStruct.accelUAV(iUAV,:,d), 'omitnan'), '-', 'LineWidth',2, 'Color',[0.83,0,0.1]);
				xlim([0 runningCorrStruct.t(end)]); % ylim([-1.1,1.1]); set(gca,'YTick',-1:0.5:1);
				box('on'); set(gca, 'FontSize',14);
				if iD==1, ylabel(['UAV #' num2str(iUAV)], 'FontSize',18); end
				if iUAV==runningCorrStruct.N, suplabel(['a_' strAx], 't', ax(:,iD), [], 'FontSize',18); end %, 'FontWeight','normal');
			end
		end
		suplabel('Time (s)', 'x', ax, [], 'FontSize',18, 'FontWeight','bold');
		suplabel('Acceleration (m/s^2)', 'y', ax, 0.03, 'FontSize',18, 'FontWeight','bold');
		legendAtBottomOfSubPlot({'Accel. (spotter''s camera)', 'Accel. (on-board IMU)'}, ax(end,end), [0.15 0.015]);
	end

	%% Plot scatterCorr: scatter a_cam vs a_UAV (for every possible combination) and fit a line
	if whatToPlot.scatterCorr
		for iD = 1:length(runningCorrStruct.dims)
			d = runningCorrStruct.dims(iD);			% For each dimension (e.g. z)
			strAx = char('X'+d-1);					% Letter representation of the dimension ('X', 'Y', 'Z')
			figs(iD, 2) = figure('Units','normalized', 'Position',[0 0 1 1]);
			ax = gobjects(runningCorrStruct.N, runningCorrStruct.M);
			for iC = 1:runningCorrStruct.M			% For each object in the cam
				iCam = runningCorrStruct.iCams(iC);
				for iUAV = 1:runningCorrStruct.N	% For each experiment
					ax(iUAV, iC) = subplot(runningCorrStruct.N,runningCorrStruct.M, iC + (iUAV-1)*runningCorrStruct.M); hold on;
					x=runningCorrStruct.accelCam(iCam,:,d); y=runningCorrStruct.accelUAV(iUAV,:,d);
					inds = union(find(isnan(x)), find(isnan(y))); x(inds)=[]; y(inds)=[];
					fitCoeffs = polyfit(x, y, 1);	% Fit a line (poly of degree 1)
					fitX = [min(x), max(x)]; fitY = polyval(fitCoeffs, fitX);
					[sigDist,p] = corr(x', y');
					hDistrib = scatter(x, y, 30); %, linspace(1,10,length(data(i).a_cam.(strAx).measured)), 'filled');
					hDistribFit = plot(fitX, fitY, '--', 'LineWidth',2, 'Color',[0.83,0,0.1]);
					text(fitX(2),fitY(2), {'\uparrow', ['y = ' num2str(fitCoeffs(1),'%.3f') 'x + ' num2str(fitCoeffs(2),'%.3f')], ['r=' num2str(sigDist,'%.3f') '; p=' num2str(p,'%.2e')]}, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',15, 'FontWeight','bold');
					box('on'); set(gca, 'FontSize',14);
					if iUAV == iCam, set(gca,'Color',correctBgndColor); end
					if iC == 1, ylabel(['a_' strAx ' IMU #' num2str(iUAV)], 'FontSize',18); end
				end
				xlabel(['a_' strAx ' cam #' num2str(iCam)], 'FontSize',18);
			end
		end
	end

	%% Plot runningLikelihoodVsWinSize: a time-domain evolution of the likelihood of each pair of cam-UAV for different likelihood computation window sizes
	if whatToPlot.runningLikelihoodVsWinSize
		figs(:, 3) = figure('Units','normalized', 'Position',[0 0 1 1]);
		ax = gobjects(runningCorrStruct.N, runningCorrStruct.M);
		for iC = 1:runningCorrStruct.M			% For each object in the cam
			iCam = runningCorrStruct.iCams(iC);
			for iUAV = 1:runningCorrStruct.N	% For each experiment
				ax(iUAV, iC) = subplot(runningCorrStruct.N,runningCorrStruct.M, iC + (iUAV-1)*runningCorrStruct.M); hold on;

				bgndNaNinds = any(isnan(runningCorrStruct.accelCam(iC,:,:)),3); startInd = 100; bgndNaNinds(1:startInd-1) = true;
				bgndT = runningCorrStruct.t(~bgndNaNinds);
				bgndY = [-5 5]; bgndY = [bgndY fliplr(bgndY)]';	% eg: [0 1 1 0]' if ylim is [0 1]

				for iW = 1:length(runningCorrStruct.frameworkWinSize)
					w = runningCorrStruct.frameworkWinSize(iW);
					fillValues = double(runningCorrStruct.assignedMatch(iUAV,startInd:end-1,iW)==iC);
					fillValues(bgndNaNinds(startInd:end-1)) = [];
					if iUAV == iCam, fprintf('Accuracy for cam #%d with a window size of %d points (%.2f s):\t%.2f%%\n', iCam, w, w*mean(diff(runningCorrStruct.t)), 100*sum(fillValues)/length(fillValues)); end

					if false % Would print the running score for each window size
						for iD = 1:length(runningCorrStruct.dims)
							d = runningCorrStruct.dims(iD);			% For each dimension (e.g. z)
							strAx = char('X'+d-1);					% Letter representation of the dimension ('X', 'Y', 'Z')
							hRunningCorr = plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningWinScore(iUAV,iC,:,iD,iW)), 'LineWidth',2);
						end
						hRunningCorr = plot(runningCorrStruct.t, squeeze(prod(runningCorrStruct.runningWinScore(iUAV,iC,:,:,iW),4)), 'LineWidth',2);
					end
					hRunningCorr = plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningPosterior(iUAV,iC,:,iW)), 'LineWidth',2);
				end
				% Let's combine contiguous patches with same fillValue together
				changeInds = find([fillValues fillValues(end)]~=[fillValues(1) fillValues]);
				newFill = ones(1,length(changeInds)+1); newFill((fillValues(1)+1):2:end) = 0;
				newBgndX = [repmat(bgndT([1, changeInds]), 2,1); repmat(bgndT([changeInds, end]), 2,1)];
				if iUAV == iCam
					if length(newFill)==3
						hBgnd = fill([newBgndX zeros(4,1)],bgndY, [newFill 0], 'FaceAlpha',0.1, 'LineStyle','None');
					else
						hBgnd = fill(newBgndX,bgndY, newFill, 'FaceAlpha',0.1, 'LineStyle','None');
					end
				else
					newBgndX(:, newFill==0) = [];	% Delete "true negatives", so we only shade in red the errors
					if ~isempty(newBgndX)
						hBgnd = fill(newBgndX,bgndY, 0, 'FaceAlpha',0.1, 'LineStyle','None');
					else
						hBgnd = [];
					end
				end
				colormap(gca, prism); set(gca, 'CLim',[0 1]);	% Prism colormap has red mapped to 0 and green to 1
				xlim([0, runningCorrStruct.t(end)]); ylim([0, 1]);
				box('on'); set(gca, 'FontSize',14, 'SortMethod','depth');	% Send background shading to bottom, so line plots are on top
				if iC==1, ylabel(['\Theta IMU #' num2str(iUAV)], 'FontSize',18); end
			end
			suplabel(['a cam obj #' num2str(iCam)], 't', ax(:,iC), [], 'FontSize',18, 'FontWeight','normal'); fprintf('\n');
		end
		suplabel('Time (s)', 'x', ax, [], 'FontSize',18);
		legendAtBottomOfSubPlot(strcat(num2str(runningCorrStruct.frameworkWinSize'), ' pt'), ax(end,end));
	end

	%% Plot runningLikelihoodFull: same as rawAccel (plot in time-domain spotter's a_cam and all a_UAV) but superimpose runningLikelihoodVsWinSize (for largest window size only) and shade based on correct assigned match
	if whatToPlot.runningLikelihoodFull
		figs(:, 4) = figure('Units','normalized', 'Position',[0 0 1 1]);
		ax = gobjects(runningCorrStruct.N, runningCorrStruct.M);
		for iC = 1:runningCorrStruct.M			% For each object in the cam
			iCam = runningCorrStruct.iCams(iC);
			for iUAV = 1:runningCorrStruct.N	% For each experiment
				ax(iUAV, iC) = subplot(runningCorrStruct.N,runningCorrStruct.M, iC + (iUAV-1)*runningCorrStruct.M); hold on;

				bgndNaNinds = any(isnan(runningCorrStruct.accelCam(iC,:,:)),3); startInd = 2; bgndNaNinds(1:startInd-1) = true;
				bgndT = runningCorrStruct.t(~bgndNaNinds);
				bgndY = [-5 5]; bgndY = [bgndY fliplr(bgndY)]';	% eg: [0 1 1 0]' if ylim is [0 1]

				for iW = min(3, length(runningCorrStruct.frameworkWinSize))
					w = runningCorrStruct.frameworkWinSize(iW);
					fillValues = double(runningCorrStruct.assignedMatch(iUAV,startInd:end-1,iW)==iC);
					fillValues(bgndNaNinds(startInd:end-1)) = [];
					if iUAV == iCam, fprintf('Accuracy for cam #%d with a window size of %d points (%.2f s):\t%.2f%%\n', iCam, w, w*mean(diff(runningCorrStruct.t)), 100*sum(fillValues)/length(fillValues)); end

					%hCam = plot(runningCorrStruct.t, runningCorrStruct.accelCam(iC,:,d), 'LineWidth',2);
					%hUAV = plot(runningCorrStruct.t, runningCorrStruct.accelUAV(iUAV,:,d), '-', 'LineWidth',2, 'Color',[0.83,0,0.1]);
					for iD = length(runningCorrStruct.dims) %%1:length(runningCorrStruct.dims)
						d = runningCorrStruct.dims(iD);			% For each dimension (e.g. z)
						strAx = char('X'+d-1);					% Letter representation of the dimension ('X', 'Y', 'Z')
						hRunningCorr = plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningWinScore(iUAV,iC,:,iD,iW)), 'LineWidth',2);
					end
					%hRunningCorr = plot(runningCorrStruct.t, squeeze(prod(runningCorrStruct.runningWinScore(iUAV,iC,:,:,iW),4)), 'LineWidth',2);
					%hRunningLikelihood = plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningLikelihood(iUAV,iC,:,iW)), 'LineWidth',2);
					hRunningPrior = plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningPosterior(iUAV,iC,:,iW)), 'LineWidth',2, 'Color',[0.83,0,0.1]);
				end
				% Let's combine contiguous patches with same fillValue together
				changeInds = find([fillValues fillValues(end)]~=[fillValues(1) fillValues]);
				newFill = ones(1,length(changeInds)+1); newFill((fillValues(1)+1):2:end) = 0;
				newBgndX = [repmat(bgndT([1, changeInds]), 2,1); repmat(bgndT([changeInds, end]), 2,1)];
				if iUAV == iC
					if length(newFill)==3
						hBgnd = fill([newBgndX zeros(4,1)],bgndY, [newFill 0], 'FaceAlpha',0.1, 'LineStyle','None');
					else
						hBgnd = fill(newBgndX,bgndY, newFill, 'FaceAlpha',0.1, 'LineStyle','None');
					end
				else
					newBgndX(:, newFill==0) = [];	% Delete "true negatives", so we only shade in red the errors
					if ~isempty(newBgndX)
						hBgnd = fill(newBgndX,bgndY, 0, 'FaceAlpha',0.1, 'LineStyle','None');
					else
						hBgnd = [];
					end
				end
				colormap(gca, prism); set(gca, 'CLim',[0 1]);	% Prism colormap has red mapped to 0 and green to 1
				xlim([0, runningCorrStruct.t(end)]); ylim([0, 1]);
				box('on'); set(gca, 'FontSize',14, 'SortMethod','depth');	% Send background shading to bottom, so line plots are on top
				if iC==1, ylabel(['a IMU #' num2str(iUAV)], 'FontSize',18); end
			end
			suplabel(['a cam obj #' num2str(iCam)], 't', ax(:,iC), [], 'FontSize',18, 'FontWeight','normal');
		end
		suplabel('Time (s)', 'x', ax, [], 'FontSize',18); fprintf('\n');
		legendAtBottomOfSubPlot({'Matching score', 'Matching confidence'}, ax(end,end));
% 		saveas(figs(1,4), 'runningLikelihoodFull', 'epsc');
% 		ch = get(figs(1,4), 'Children');
% 		for i = 1:length(ch)
% 			if isa(ch(i),'matlab.graphics.axis.Axes')
% 				ch2 = get(ch(i), 'Children');
% 				for j = 1:length(ch2)
% 					set(ch2(j), 'Visible','off');
% 				end
% 			end
% 		end
% 		saveas(figs(1,4), 'runningLikelihoodFull_empty', 'epsc');
	end
end
