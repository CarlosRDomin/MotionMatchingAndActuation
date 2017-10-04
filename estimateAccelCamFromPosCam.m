%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applies the Savitzky-Golay derivative filter to a time-window of the past frameworkWinSize points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [accelCam] = estimateAccelCamFromPosCam(posUAVcam, accelCam, currTind, derivFiltOrder, derivFiltHalfWinSize, frameworkWinSize, fps)
	% First and last derivFiltHalfWinSize points are incorrect (window effects). Which means 2 things:
	%  - The time window should be shifted back by derivFiltHalfWinSize points
	%  - For the point at currTind-derivFiltHalfWinSize to be correct, we need to start computing the filter at currTind-2*derivFiltHalfWinSize
	tDerivFilterInds = currTind+(-2*derivFiltHalfWinSize:frameworkWinSize);	tDerivFilterInds(tDerivFilterInds<1) = []; % Make the first time-windows don't yield "index out of bounds" errors
	
	% Save the results of the filter applied at tDerivFilterInds, we will then index it and use from the derivFiltHalfWinSize+1'th point on
	filteredAccelCam = derivFilter(posUAVcam(:,tDerivFilterInds,:), 2, 1/fps, derivFiltOrder, 2*derivFiltHalfWinSize+1);
	%if currTind > 1
		validDerivFilterInds = derivFiltHalfWinSize+1 : length(tDerivFilterInds)-derivFiltHalfWinSize; % Remember to avoid the derivFiltHalfWinSize points on either end of the window
	%else
	%	validDerivFilterInds = 1:length(tDerivFilterInds);
	%end
	
	% Finally, save the results (in the valid indices)
	accelCam(:,tDerivFilterInds(validDerivFilterInds),:) = filteredAccelCam(:,validDerivFilterInds,:);
end
