function generateDronesInRoomVideo(videoName, posUAVgt, roomDimensions, spotterCam)
	if isempty(videoName), videoName = 'videoSim.mp4'; end
	tempFrameName = 'frame.jpg';

	% Initialize the video object
	v = VideoWriter(videoName,'MPEG-4');
	open(v);
	
	% For the first frame, need to create the figure and the handler
	figure('Units','normalized', 'Position',[0.3 0.4 0.4 0.25]);
	h=plotDronesInRoom(posUAVgt(:,1,:), roomDimensions, spotterCam);
	saveas(gcf, tempFrameName);
	writeVideo(v,imread(tempFrameName));
	
	% For all other frames, simply update XData, YData and ZData
	for k=2:5:size(posUAVgt,2)
		dispImproved(sprintf('\nProcessing video frame %4d out of %4d (%6.2f%%)', k, size(posUAVgt,2), 100*k/size(posUAVgt,2)));
		%pause(1/(2*pointsPerIter)); % Play at 2x "real-time" 
		set(h, 'XData',posUAVgt(:,k,1), 'YData',posUAVgt(:,k,2), 'ZData',posUAVgt(:,k,3));
		saveas(gcf, tempFrameName);
		writeVideo(v,imread(tempFrameName));
	end
	
	% Done! :)
	close(v);
	delete(tempFrameName);
	dispImproved(sprintf('\nDone saving simulation video!\n'), 'keepthis');
end
