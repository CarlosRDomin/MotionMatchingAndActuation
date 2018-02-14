function plotAccuracyAndSurvivorVsN(experimentsStruct, typeOfMotion, numWindowsPerSecond, tEnd, intervalErrorBar, lWidth, lWidthErr, fSizeLabels, fSizeAxes)
	if nargin<3 || isempty(numWindowsPerSecond)
		numWindowsPerSecond = 1;
	end
	if nargin<4 || isempty(tEnd)
		tEnd = 20;
	end
	if nargin<5 || isempty(intervalErrorBar)
		intervalErrorBar = 2;
	end
	if nargin<6 || isempty(lWidth)
		lWidth = 2;
	end
	if nargin<7 || isempty(lWidthErr)
		lWidthErr = 1;
	end
	if nargin<8 || isempty(fSizeLabels)
		fSizeLabels = 17;
	end
	if nargin<9 || isempty(fSizeAxes)
		fSizeAxes = 13.5;
	end
	
	figure('Units','pixels', 'Position',[200 200, 560 375]);
	ax = gobjects(2,1); h = gobjects(2,numel(experimentsStruct.(typeOfMotion).N));

	% Compute MOTA and numSurvivors stats (mean, std)
	[MOTAmean, MOTAstd, numSurvivorsMean, numSurvivorsStd] = computeMotaAndSurvivalStats(experimentsStruct, typeOfMotion, intervalErrorBar);
	t = (0:size(numSurvivorsMean,2)-1) / numWindowsPerSecond;

	for Nind = 1:numel(experimentsStruct.(typeOfMotion).N)
% 		% Plot MOTA
% 		ax(1) = subplot(2,1,1); hold on;
% 		h(1,Nind) = errorbar(t, MOTAmean(Nind,:), MOTAstd(Nind,:), 'LineWidth',lWidthErr);
% 		plot(t, MOTAmean(Nind,:), 'LineWidth',lWidth, 'Color',get(h(1,Nind),'Color'));
% 		xlim([0,tEnd]); ylim([0 100]); ylabel('Accuracy (%)', 'FontSize',fSizeLabels+1);
% 
% 		% Plot survival rate
% 		ax(2) = subplot(2,1,2); hold on;
% 		h(2,Nind) = errorbar(t, numSurvivorsMean(Nind,:), numSurvivorsStd(Nind,:), 'LineWidth',lWidthErr);
% 		plot(t, numSurvivorsMean(Nind,:), 'LineWidth',lWidth, 'Color',get(h(2,Nind),'Color'));
% 		xlim([0,tEnd]); ylim([-0.05 1.05]); ylabel('Survival rate (%)', 'FontSize',fSizeLabels+1);
		[h(:,Nind), ax] = plotAccuracyAndSurvivor(t, MOTAmean, MOTAstd, numSurvivorsMean, numSurvivorsStd, Nind, tEnd, lWidth, lWidthErr, fSizeLabels);
	end

	suplabel('Time (s)', 'x', ax(:), -0.01, 'FontSize',fSizeLabels);
	set(ax(:), 'Box','on', 'FontSize',fSizeAxes);
	l=legend(h(1,:), cellstr(strcat('N=',num2str((experimentsStruct.(typeOfMotion).N)'))), 'Orientation','horizontal', 'FontSize',fSizeAxes); set(l, 'Units','normalized', 'Position',[0.17 0.494 0.7 0]);
end
