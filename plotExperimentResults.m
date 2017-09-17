%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates different plots (specified by the struct whatToPlot) from a set of real experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [runningCorrStruct] = plotExperimentResults(runningCorrStruct, whatToPlot, data)
	%% Load parameters / default values
	if nargin<1 || isempty(runningCorrStruct)
		try
			aux = load('runningCorrStruct.mat');
			runningCorrStruct = aux.runningCorrStruct;
			runningCorrStruct.runningCorrWinSizes([1:2 4:end]) = [];
			runningCorrStruct.runningWinScore(:,:,:,:,[1:2 4:end]) = [];
			runningCorrStruct.assignedMatch(:,:,:,[1:2 4:end]) = [];
			runningCorrStruct.runningPrior(:,:,:,:,[1:2 4:end]) = [];
			runningCorrStruct.runningLikelihood(:,:,:,:,[1:2 4:end]) = [];
			runningCorrStruct.dims = 3;
			error;	% Force to go into loading experiment data from scratch
		catch
			if nargin<3 || isempty(data)	% Only load data if it hasn't been loaded yet (saves time)
				%data = loadExperimentData(struct('datetime',{'2017-02-19 17-14-15','2017-02-19 17-56-48','2017-02-19 17-59-23','2017-02-19 18-01-29','2017-02-19 18-22-23','2017-02-19 18-24-25'}, 'ch','75'));
				data = loadExperimentData(struct('datetime',{'2017-02-19 17-56-48','2017-02-19 17-59-23','2017-02-19 18-01-29','2017-02-19 18-22-23','2017-02-19 18-24-25'}, 'ch','75'));
			end
			runningCorrStruct = runMatchingFrameworkOnGivenData(data, [1 2 3], [1 2 4], [5:5:15 30:30:300]);	%(data, iCam=1:length(data), runningCorrWinSizes, dims=1:3, cropToMinT=true)
			save('runningCorrStruct.mat','runningCorrStruct');
		end
	end
	if nargin<2 || isempty(whatToPlot)
		whatToPlot = struct('rawAccel',false, 'scatterCorr',false, 'runningLikelihoodVsWinSize',false, 'runningLikelihoodFull',true);
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
					hCam = plot(runningCorrStruct.t, runningCorrStruct.yCam(runningCorrStruct.iCams==iUAV,:,d), 'LineWidth',2);
				end
				hUAV = plot(runningCorrStruct.t, runningCorrStruct.yUAV(iUAV,:,d)-mean(runningCorrStruct.yUAV(iUAV,:,d), 'omitnan'), '-', 'LineWidth',2, 'Color',[0.83,0,0.1]);
				xlim([0 runningCorrStruct.t(end)]); % ylim([-1.1,1.1]); set(gca,'YTick',-1:0.5:1);
				box('on'); set(gca, 'FontSize',14);
				if iD==1, ylabel(['UAV #' num2str(iUAV)], 'FontSize',18); end
				if iUAV==runningCorrStruct.N, suplabel(['a_' strAx], 't', ax(:,iD), [], 'FontSize',18); end %, 'FontWeight','normal');
			end
		end
		suplabel('Time (s)', 'x', ax, [], 'FontSize',18, 'FontWeight','bold');
		suplabel('Acceleration (m/s^2)', 'y', ax, 0.03, 'FontSize',18, 'FontWeight','bold');
		legendAtBottomOfSubPlot({['a_' strAx ' (spotter''s camera)'], ['a_' strAx ' (on-board IMU)']}, ax(end,end), [0.15 0.015]);
		% f=gcf; f.Children = flip(f.Children);
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
					x=runningCorrStruct.yCam(iCam,:,d); y=runningCorrStruct.yUAV(iUAV,:,d);
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

				bgndNaNinds = any(isnan(runningCorrStruct.yCam(iC,:,:)),3); startInd = 100; bgndNaNinds(1:startInd-1) = true;
				bgndT = runningCorrStruct.t(~bgndNaNinds);
				bgndY = [-5 5]; bgndY = [bgndY fliplr(bgndY)]';	% eg: [0 1 1 0]' if ylim is [0 1]

				for iW = 1:length(runningCorrStruct.runningCorrWinSizes)
					w = runningCorrStruct.runningCorrWinSizes(iW);
					fillValues = double(runningCorrStruct.assignedMatch(iUAV,startInd:end-1,iW)==iC);
					fillValues(bgndNaNinds(startInd:end-1)) = [];
					if iUAV == iCam, fprintf('Accuracy for cam #%d with a window size of %d points (%.2f s):\t%.2f%%\n', iCam, w, w*mean(diff(runningCorrStruct.t)), 100*sum(fillValues)/length(fillValues)); end

% 					for iD = 1:length(runningCorrStruct.dims)
% 						d = runningCorrStruct.dims(iD);			% For each dimension (e.g. z)
% 						strAx = char('X'+d-1);					% Letter representation of the dimension ('X', 'Y', 'Z')
% 						hRunningCorr = plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningWinScore(iUAV,iC,:,iD,iW)), 'LineWidth',2);
% 					end
% 					hRunningCorr = plot(runningCorrStruct.t, squeeze(prod(runningCorrStruct.runningWinScore(iUAV,iC,:,:,iW),4)), 'LineWidth',2);
					hRunningCorr = plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningPrior(iUAV,iC,2:end,iW)), 'LineWidth',2);
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
				%if iC==1, ylabel(['\Theta_' strAx ' IMU #' num2str(iUAV)], 'FontSize',18); end
				if iC==1, ylabel(['\Theta IMU #' num2str(iUAV)], 'FontSize',18); end
			end
			%suplabel(['a_' strAx ' cam obj #' num2str(iCam)], 't', ax(:,iC), [], 'FontSize',18, 'FontWeight','normal'); fprintf('\n');
			suplabel(['a cam obj #' num2str(iCam)], 't', ax(:,iC), [], 'FontSize',18, 'FontWeight','normal'); fprintf('\n');
		end
		suplabel('Time (s)', 'x', ax, [], 'FontSize',18);
		legendAtBottomOfSubPlot(strcat(num2str(runningCorrStruct.runningCorrWinSizes'), ' pt'), ax(end,end));
	end
	%% Plot runningLikelihoodFull: same as rawAccel (plot in time-domain spotter's a_cam and all a_UAV) but superimpose runningLikelihoodVsWinSize (for largest window size only) and shade based on correct assigned match
	if whatToPlot.runningLikelihoodFull
		figs(:, 4) = figure('Units','normalized', 'Position',[0 0 1 1]);
		ax = gobjects(runningCorrStruct.N, runningCorrStruct.M);
		for iC = 1:runningCorrStruct.M			% For each object in the cam
			iCam = runningCorrStruct.iCams(iC);
			for iUAV = 1:runningCorrStruct.N	% For each experiment
				ax(iUAV, iC) = subplot(runningCorrStruct.N,runningCorrStruct.M, iC + (iUAV-1)*runningCorrStruct.M); hold on;

				bgndNaNinds = any(isnan(runningCorrStruct.yCam(iC,:,:)),3); startInd = 2; bgndNaNinds(1:startInd-1) = true;
				bgndT = runningCorrStruct.t(~bgndNaNinds);
				bgndY = [-5 5]; bgndY = [bgndY fliplr(bgndY)]';	% eg: [0 1 1 0]' if ylim is [0 1]

				for iW = min(3, length(runningCorrStruct.runningCorrWinSizes))
					w = runningCorrStruct.runningCorrWinSizes(iW);
					fillValues = double(runningCorrStruct.assignedMatch(iUAV,startInd:end-1,iW)==iC);
					fillValues(bgndNaNinds(startInd:end-1)) = [];
					if iUAV == iCam, fprintf('Accuracy for cam #%d with a window size of %d points (%.2f s):\t%.2f%%\n', iCam, w, w*mean(diff(runningCorrStruct.t)), 100*sum(fillValues)/length(fillValues)); end

					%hCam = plot(runningCorrStruct.t, runningCorrStruct.yCam(iC,:,d), 'LineWidth',2);
					%hUAV = plot(runningCorrStruct.t, runningCorrStruct.yUAV(iUAV,:,d), '-', 'LineWidth',2, 'Color',[0.83,0,0.1]);
					for iD = 3 %%1:length(runningCorrStruct.dims)
						d = runningCorrStruct.dims(iD);			% For each dimension (e.g. z)
						strAx = char('X'+d-1);					% Letter representation of the dimension ('X', 'Y', 'Z')
						hRunningCorr = plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningWinScore(iUAV,iC,:,iD,iW)), 'LineWidth',2);
					end
					%%hRunningCorr = plot(runningCorrStruct.t, squeeze(prod(runningCorrStruct.runningWinScore(iUAV,iC,:,:,iW),4)), 'LineWidth',2);
					%%hRunningLikelihood = plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningLikelihood(iUAV,iC,:,iW)), 'LineWidth',2);
					hRunningPrior = plot(runningCorrStruct.t, squeeze(runningCorrStruct.runningPrior(iUAV,iC,2:end,iW)), 'LineWidth',2, 'Color',[0.83,0,0.1]);
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
				%if iC==1, ylabel(['a_' strAx ' IMU #' num2str(iUAV)], 'FontSize',18); end
				if iC==1, ylabel(['a IMU #' num2str(iUAV)], 'FontSize',18); end
				% if iC==1, ylabel(['a_' strAx ' IMU #' char('A'+iUAV-1)], 'FontSize',18); end
			end
			%suplabel(['a_' strAx ' cam obj #' num2str(iCam)], 't', ax(:,iC), [], 'FontSize',18, 'FontWeight','normal');
			suplabel(['a cam obj #' num2str(iCam)], 't', ax(:,iC), [], 'FontSize',18, 'FontWeight','normal');
		end
		suplabel('Time (s)', 'x', ax, [], 'FontSize',18); fprintf('\n');
		%legendAtBottomOfSubPlot({['a_' strAx ' (spotter''s camera #' num2str(iCam) ')'], ['a_' strAx ' (on-board IMU)'], 'Running score', 'Running likelihood', 'Running confidence'});
		%%legendAtBottomOfSubPlot({'S_X', 'S_Y', 'S_Z', 'Matching score', 'Matching likelihood', 'Matching confidence'}, ax(end,end));
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
