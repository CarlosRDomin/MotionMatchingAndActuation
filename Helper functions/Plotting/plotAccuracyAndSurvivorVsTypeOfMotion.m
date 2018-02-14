function plotAccuracyAndSurvivorVsTypeOfMotion(experimentsStruct, typeOfMotionCell, legendCell, N, numWindowsPerSecond, tEnd, intervalErrorBar, lWidth, lWidthErr, fSizeLabels, fSizeAxes)
	if nargin<5 || isempty(numWindowsPerSecond)
		numWindowsPerSecond = 1;
	end
	if nargin<6 || isempty(tEnd)
		tEnd = 20;
	end
	if nargin<7 || isempty(intervalErrorBar)
		intervalErrorBar = 2;
	end
	if nargin<8 || isempty(lWidth)
		lWidth = 2;
	end
	if nargin<9 || isempty(lWidthErr)
		lWidthErr = 1;
	end
	if nargin<10 || isempty(fSizeLabels)
		fSizeLabels = 17;
	end
	if nargin<11 || isempty(fSizeAxes)
		fSizeAxes = 13.5;
	end
	
	figure('Units','pixels', 'Position',[200 200, 560 375]);
	h = gobjects(2,numel(typeOfMotionCell));
	for typeOfMotionInd = 1:numel(typeOfMotionCell)
		typeOfMotion = typeOfMotionCell{typeOfMotionInd};
		[~,Nind] = find(experimentsStruct.(typeOfMotion).N==N);
			
		% Compute accuracy and numSurvivors stats (mean, std)
		[MOTAmean, MOTAstd, numSurvivorsMean, numSurvivorsStd] = computeMotaAndSurvivalStats(experimentsStruct, typeOfMotion, intervalErrorBar);
		t = (0:size(numSurvivorsMean,2)-1) / numWindowsPerSecond;

		[h(:,typeOfMotionInd), ax] = plotAccuracyAndSurvivor(t, MOTAmean, MOTAstd, numSurvivorsMean, numSurvivorsStd, Nind, tEnd, lWidth, lWidthErr, fSizeLabels);
	end
	suplabel('Time (s)', 'x', ax(:), -0.01, 'FontSize',fSizeLabels);
	set(ax(:), 'Box','on', 'FontSize',fSizeAxes);
	l=legend(h(1,:), legendCell, 'Orientation','horizontal', 'FontSize',fSizeAxes, 'Interpreter','latex'); set(l, 'Units','normalized', 'Position',[0.17 0.494 0.7 0]);
end
