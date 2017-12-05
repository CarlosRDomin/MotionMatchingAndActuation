function [h, ax] = plotDronesInRoom(posUAVgt, roomDimensions, spotterCamInfo, scatterPtSize)
	if nargin<3, spotterCamInfo = []; end
	if nargin<4, scatterPtSize = []; end
	ax = zeros(1,3); h = zeros(1,3);
	spotterCamInfoCopy = spotterCamInfo;
	spotterCamInfoCopy.alpha = 0.05;
	
	ax(1) = subplot(1,3,1); h(1) = plotDronesInGCA(posUAVgt, roomDimensions, spotterCamInfo, scatterPtSize); view([90, 90]); title('Top view');
	ax(2) = subplot(1,3,2); h(2) = plotDronesInGCA(posUAVgt, roomDimensions, spotterCamInfoCopy, scatterPtSize); view([90, 0]); title('Front view');
	ax(3) = subplot(1,3,3); h(3) = plotDronesInGCA(posUAVgt, roomDimensions, spotterCamInfo, scatterPtSize); view([0, 0]); title('Side view');
	
	% Reduce margins and set subplot width proportional to their horizontal limits (xlim() or ylim() depending on view) -> Same scale
	marginH = [0.05, 0.03]; marginV = [0.125, 0.075]; padding = 0.05;
	totalHorizAxes = sum(diff([ylim(); ylim(); xlim()]'));
	widthAxes = [diff(ylim()); diff(ylim()); diff(xlim())]/totalHorizAxes*(1-sum(marginH)-2*padding);
	posStartAxes = marginH(1) + [0; cumsum(widthAxes(1:end-1))] + (0:2)'*padding;
	newPosAxes = [posStartAxes repmat(marginV(1), 3,1) widthAxes repmat(1-sum(marginV), 3,1)];
	for i=1:3, set(ax(i), 'Position', newPosAxes(i,:)); end
end

function h = plotDronesInGCA(posUAVgt, roomDimensions, spotterCamInfo, scatterPtSize)
	maxLims = roomDimensions;
	if nargin<3 || isempty(spotterCamInfo)
		spotterCamInfo = [];
	else
		maxLims(1) = spotterCamInfo.pos(1)+0.5;
		if ~isfield(spotterCamInfo, 'size'), spotterCamInfo.size = 0.25.*[1, 3/4, 1]; end
		if ~isfield(spotterCamInfo, 'alpha'), spotterCamInfo.alpha = 0.15; end
	end
	if nargin<4 || isempty(scatterPtSize), scatterPtSize = 80; end

	cla; hold on; colormap('jet');
	h = scatter3(posUAVgt(:,1,1), posUAVgt(:,1,2), posUAVgt(:,1,3), scatterPtSize, (1:size(posUAVgt,1))', 'filled');
	if ~isempty(spotterCamInfo), plotCamera('Location',spotterCamInfo.pos, 'Orientation',spotterCamInfo.orient, 'Size',spotterCamInfo.size, 'Opacity',spotterCamInfo.alpha); end
	plotCube(zeros(1,3), roomDimensions, 'FaceColor','blue', 'EdgeColor','blue', 'FaceAlpha',.1, 'EdgeAlpha',.5);
	grid('on'); axis equal; xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
	lims = [zeros(1,3); maxLims]'; xlim(lims(1,:)); xticks(lims(1,1):lims(1,2)); ylim(lims(2,:)); yticks(floor(lims(2,1)):lims(2,2)); zlim(lims(3,:)); zticks(lims(3,1):lims(3,2));
end
