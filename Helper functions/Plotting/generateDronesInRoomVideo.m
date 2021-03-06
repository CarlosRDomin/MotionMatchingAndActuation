function generateDronesInRoomVideo(videoName, posUAVgtOrExperimentInfo, roomDimensions, spotterCam, saveVideo)
	if isempty(videoName), videoName = 'videoSim.mp4'; end
	if isstruct(posUAVgtOrExperimentInfo) % Allow to visualize data directly from a simulation *.mat
		if isempty(videoName), [~,videoName] = fileparts(posUAVgtOrExperimentInfo.experimentMatFileName); end % Default output filename based on simulation filename
		[expData, experimentInd] = loadSimulationExperimentData(posUAVgtOrExperimentInfo); % Load simulation data
		posUAVgt = expData.variableStruct(experimentInd).posUAVgt; % Populate fields from the data
		roomDimensions = expData.paramStruct.roomDimensions;
		spotterCam = expData.paramStruct.spotterCam;
		spotterCam.orient = [0 1 0; 0 0 -1; -1 0 0];
		spotterCam.pos = [roomDimensions(1) + (roomDimensions(1)/2) / tan(spotterCam.FOV/2), roomDimensions(2)/2, roomDimensions(3)/2];
	else % Otherwise, posUAVgtOrExperimentInfo is in fact the posUAVgt to visualize
		posUAVgt = posUAVgtOrExperimentInfo;
	end
	if nargin<5 || isempty(saveVideo), saveVideo = true; end
	tempFrameName = 'frame.jpg';

	% Initialize the video object
	if saveVideo
		v = VideoWriter(videoName,'MPEG-4');
		v.FrameRate = spotterCam.fps;
		open(v);
	end
	
	% For the first frame, need to create the figure and the handler
	figure('Units','normalized', 'Position',[0.3 0.4 0.4 0.25]);
	[h,ax] = plotDronesInRoom(posUAVgt(:,1,:), roomDimensions, spotterCam);
	[~,hTitle]=suplabel(['Time: ' num2str(0, '%.2f')], 't', ax(:), 0);
	
	% For all other frames, simply update XData, YData and ZData
	for k=1:size(posUAVgt,2)
		dispImproved(sprintf('Processing video frame %4d out of %4d (%6.2f%%)', k, size(posUAVgt,2), 100*k/size(posUAVgt,2)));
		pause(1/spotterCam.fps);
		set(h, 'XData',posUAVgt(:,k,1), 'YData',posUAVgt(:,k,2), 'ZData',posUAVgt(:,k,3));
		set(hTitle, 'String',['Time: ' num2str(k/spotterCam.fps, '%.2f')]);
		if saveVideo
			saveas(gcf, tempFrameName);
			writeVideo(v,imread(tempFrameName));
		end
	end
	
	% Done! :)
	if saveVideo
		delete(tempFrameName);
		close(v);
	end
	dispImproved(sprintf('Done saving simulation video!\n'), 'keepthis');
end
