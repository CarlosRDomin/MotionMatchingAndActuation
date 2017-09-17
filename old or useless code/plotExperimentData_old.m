close all;
includeScripts;

if exist('data','var') ~= 1	% Only load data if it hasn't been loaded yet (saves time)
	data = loadExperimentData(struct('datetime',{'2017-02-19 17-14-15','2017-02-19 17-56-48','2017-02-19 17-59-23','2017-02-19 18-01-29','2017-02-19 18-22-23','2017-02-19 18-24-25'}, 'ch','75'));
end
% data = loadExperimentData(struct('datetime',{'2017-02-19 17-56-48'}, 'ch','75'));
whatToPlot = struct('rawAccel',true, 'rawGyro',false, 'fullCross',true, 'dataDistrib',true);
figs = zeros(3, 1); % 3 stands for {X,Y,Z}
figSize = [560, 250];

for i = 1:length(data)  % For each experiment
    for d = 1:3         % For each dimension (x, y, z)
        strAx = char('X'+d-1);  % Letter representation of the dimension ('X', 'Y', 'Z')
        
		if whatToPlot.rawAccel
			figs(d) = figure('Units','pixels', 'Position',[(d-1)*figSize(1),200, figSize]); hold on;
			hCam = plot(data(i).a_cam.(strAx).t, data(i).a_cam.(strAx).measured, 'LineWidth',2);
			hUAV = plot(data(i).a_UAV.(strAx).t, data(i).a_UAV.(strAx).measured, '-', 'LineWidth',2, 'Color',[0.83,0,0.1]);
			% set(hCam, 'XData',get(hCam, 'XData') - xl(1)); set(hUAV, 'XData',get(hUAV, 'XData') - xl(1)); xlim(xl-xl(1)); % Change xlim, then update x-axis values so time starts at t=0
			% ylim([-1.1,1.1]); set(gca,'YTick',-1:0.5:1);
			box('on'); set(gca, 'FontSize',14);
			xlabel('Time (s)', 'FontSize',18);
			ylabel('Acceleration (m/s^2)', 'FontSize',18);
			l=legend(['a_' strAx ' (spotter''s camera)'], ['a_' strAx ' (on-board IMU)'], 'Location','SouthWest'); set(l, 'FontSize',13);
			%%saveas(figs(i), ['a_' strAx], 'epsc');
		end

		if whatToPlot.rawGyro
			figure('Units','pixels', 'Position',[(d-1)*figSize(1),200+2*figSize(2), figSize]);
			hGyro = plot(data(i).gyro_UAV.(strAx).t, data(i).gyro_UAV.(strAx).measured, '-', 'LineWidth',2, 'Color',[0.83,0,0.1]);
			box('on'); set(gca, 'FontSize',14); xlabel('Time (s)', 'FontSize',18); ylabel('Gyro', 'FontSize',18); l=legend(['gyro_' strAx ' (on-board IMU)'], 'Location','SouthWest'); set(l, 'FontSize',13);
		end

		if whatToPlot.fullCross
			figure('Units','pixels', 'Position',[(d-1)*figSize(1),200+figSize(2), figSize]);
			T_s = mean(diff(data(i).a_cam.(strAx).t));
			[crossCorr, lags] = xcorr(data(i).a_cam.(strAx).measured, data(i).a_UAV.(strAx).measured, round(2/T_s));  % Max lag 2sec
			hCross = stem(lags*T_s, crossCorr, '-', 'LineWidth',2);
			box('on'); set(gca, 'FontSize',14); xlabel('Lag (s)', 'FontSize',18); ylabel('Cross-correlation', 'FontSize',18);
		end

		if whatToPlot.dataDistrib
			figure('Units','pixels', 'Position',[(d-1)*figSize(1),200+3*figSize(2), figSize]); hold on;
			fitCoeffs = polyfit(data(i).a_cam.(strAx).measured, data(i).a_UAV.(strAx).measured, 1);	% Fit a line
			fitX = [min(data(i).a_cam.(strAx).measured), max(data(i).a_cam.(strAx).measured)]; fitY = polyval(fitCoeffs, fitX);
			[R,p] = corr(data(i).a_cam.(strAx).measured, data(i).a_UAV.(strAx).measured);
			hDistrib = scatter(data(i).a_cam.(strAx).measured, data(i).a_UAV.(strAx).measured, 30); %, linspace(1,10,length(data(i).a_cam.(strAx).measured)), 'filled');
			hDistribFit = plot(fitX, fitY, '--', 'LineWidth',2, 'Color',[0.83,0,0.1]); text(fitX(2),fitY(2), {'\uparrow', ['y = ' num2str(fitCoeffs(1),'%.3f') 'x + ' num2str(fitCoeffs(2),'%.3f')], ['r=' num2str(R,'%.3f') '; p=' num2str(p,'%.2e')]}, 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',15, 'FontWeight','bold');
			box('on'); set(gca, 'FontSize',14); xlabel('Cam Acceleration (m/s^2)', 'FontSize',18); ylabel('UAV Acceleration (m/s^2)', 'FontSize',18);
		end
    end
end