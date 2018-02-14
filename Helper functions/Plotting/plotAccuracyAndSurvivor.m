function [h, ax] = plotAccuracyAndSurvivor(t, MOTAmean, MOTAstd, numSurvivorsMean, numSurvivorsStd, Nind, tEnd, lWidth, lWidthErr, fSizeLabels)
	if nargin<8 || isempty(lWidth)
		lWidth = 2;
	end
	if nargin<9 || isempty(lWidthErr)
		lWidthErr = 1;
	end
	if nargin<10 || isempty(fSizeLabels)
		fSizeLabels = 17;
	end
	h = gobjects(2,1); ax = gobjects(2,1);
	
	% Plot accuracy
	ax(1) = subplot(2,1,1); hold on;
	h(1) = errorbar(t, MOTAmean(Nind,:), MOTAstd(Nind,:), 'LineWidth',lWidthErr);
	plot(t, MOTAmean(Nind,:), 'LineWidth',lWidth, 'Color',get(h(1),'Color'));
	xlim([0,tEnd]); ylim([0 100]); ylabel('Accuracy (%)', 'FontSize',fSizeLabels+1);

	% Plot survival rate
	ax(2) = subplot(2,1,2); hold on;
	h(2) = errorbar(t, numSurvivorsMean(Nind,:), numSurvivorsStd(Nind,:), 'LineWidth',lWidthErr);
	plot(t, numSurvivorsMean(Nind,:), 'LineWidth',lWidth, 'Color',get(h(2),'Color'));
	xlim([0,tEnd]); ylim([0 100]); ylabel('Survival rate (%)', 'FontSize',fSizeLabels+1);
end
