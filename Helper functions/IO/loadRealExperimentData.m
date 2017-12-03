%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads (real, ie. from CF drones) experiment data from .npz files into a struct which contains the accel measured from the camera, accel and gyro measured from the drone's IMU...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = loadRealExperimentData(logInfo, logFolder, derivPolyOrder, derivWinSize, movingAvgFiltWinSize)
	if nargin<2 || isempty(logFolder), logFolder = [fileparts(mfilename('fullpath')) '/../../data/Real']; end
	if nargin<3 || isempty(derivPolyOrder), derivPolyOrder = 2; end
	if nargin<4 || isempty(derivWinSize), derivWinSize = 1 + 2*10; end
	if nargin<5 || isempty(movingAvgFiltWinSize), movingAvgFiltWinSize = 30; end
    magnitudes = {'X'; 'Y'; 'Z'};
    
    data = repmat(cell2struct([repmat({cell2struct(cell(3,1), magnitudes)}, 6,1); cell(2,1)], {'a_cam','a_UAV','gyro_UAV','a_cam_orig','a_UAV_orig','gyro_UAV_orig', 'tInterv', 'tCropInds'}), length({logInfo.datetime}), 1);
    for j = 1:length({logInfo.datetime})
		drone_id_filename = [logFolder '/' logInfo(j).datetime '/log_droneId_' logInfo(j).ch '_' strrep(logInfo(j).datetime, ' ','_') '.npz'];
		if exist(drone_id_filename, 'file'), drone_id_orig = readNPZ(drone_id_filename); else, drone_id_orig = []; end

        for i = 1:size(magnitudes,1)
			% Read logs from files
			pos_cam_orig = readNPZ([logFolder '/' logInfo(j).datetime '/log_p' magnitudes{i} '_cam_' logInfo(j).ch '_' strrep(logInfo(j).datetime, ' ','_') '.npz']);
			vel_cam_orig = readNPZ([logFolder '/' logInfo(j).datetime '/log_v' magnitudes{i} '_cam_' logInfo(j).ch '_' strrep(logInfo(j).datetime, ' ','_') '.npz']);
            accel_cam_orig = readNPZ([logFolder '/' logInfo(j).datetime '/log_a' magnitudes{i} '_cam_' logInfo(j).ch '_' strrep(logInfo(j).datetime, ' ','_') '.npz']);
            accel_UAV_orig = readNPZ([logFolder '/' logInfo(j).datetime '/log_a' magnitudes{i} '_world_' logInfo(j).ch '_' strrep(logInfo(j).datetime, ' ','_') '.npz']); accel_UAV_orig.measured = 9.81.*(accel_UAV_orig.measured - accel_UAV_orig.measured(1));
            gyro_UAV_orig = readNPZ([logFolder '/' logInfo(j).datetime '/log_gyro' magnitudes{i} '_' logInfo(j).ch '_' strrep(logInfo(j).datetime, ' ','_') '.npz']);
            accel_UAV_dsamp = struct('tFloat',accel_cam_orig.tFloat, 'measured',interp1(accel_UAV_orig.tFloat, movingAvgFilter(movingAvgFiltWinSize, accel_UAV_orig.measured), accel_cam_orig.tFloat)); % accel_UAV_dsamp.measured(isnan(accel_UAV_dsamp.measured))=0;
            gyro_UAV_dsamp = struct('tFloat',accel_cam_orig.tFloat, 'measured',interp1(gyro_UAV_orig.tFloat, movingAvgFilter(movingAvgFiltWinSize, gyro_UAV_orig.measured), accel_cam_orig.tFloat)); % gyro_UAV_dsamp.measured(isnan(gyro_UAV_dsamp.measured))=0;
			
			% Apply the derivative filter to convert pos_cam_orig.measured to vel and accel
			vel_cam_orig.measured = derivFilter(reshape(pos_cam_orig.measured, 1,[]), 1, mean(diff(pos_cam_orig.tFloat)), derivPolyOrder, derivWinSize);
			accel_cam_orig.measured = derivFilter(reshape(pos_cam_orig.measured, 1,[]), 2, mean(diff(pos_cam_orig.tFloat)), derivPolyOrder, derivWinSize);

			% Crop the variables according to accel_cam_orig.tFloat
            xl = [accel_cam_orig.tFloat(1), accel_UAV_orig.tFloat(end)];
            %xl = [accel_cam_orig.tFloat((derivWinSize+1)/2), accel_UAV_orig.tFloat(end-(3*derivWinSize-1)/2)];
            tCropInds = (accel_cam_orig.tFloat>=xl(1) & accel_cam_orig.tFloat<=xl(2));
            tCropIndsUAV = (accel_UAV_orig.tFloat>=xl(1) & accel_UAV_orig.tFloat<=xl(2));
			if ~isempty(drone_id_orig), drone_id_crop = struct('tFloat',drone_id_orig.tFloat(tCropInds), 't',drone_id_orig.tFloat(tCropInds)-drone_id_orig.tFloat(find(tCropInds,1)) , 'measured',drone_id_orig.measured(tCropInds)); end
			pos_cam_crop = struct('tFloat',pos_cam_orig.tFloat(tCropInds), 't',pos_cam_orig.tFloat(tCropInds)-pos_cam_orig.tFloat(find(tCropInds,1)) , 'measured',pos_cam_orig.measured(tCropInds));
			vel_cam_crop = struct('tFloat',vel_cam_orig.tFloat(tCropInds), 't',vel_cam_orig.tFloat(tCropInds)-vel_cam_orig.tFloat(find(tCropInds,1)) , 'measured',vel_cam_orig.measured(tCropInds));
			accel_cam_crop = struct('tFloat',accel_cam_orig.tFloat(tCropInds), 't',accel_cam_orig.tFloat(tCropInds)-accel_cam_orig.tFloat(find(tCropInds,1)) , 'measured',accel_cam_orig.measured(tCropInds));
            accel_UAV_crop = struct('tFloat',accel_UAV_orig.tFloat(tCropIndsUAV), 't',accel_UAV_orig.tFloat(tCropIndsUAV)-accel_UAV_orig.tFloat(find(tCropIndsUAV,1)), 'measured',accel_UAV_orig.measured(tCropIndsUAV));
			accel_UAV_dsamp_crop = struct('tFloat',accel_UAV_dsamp.tFloat(tCropInds), 't',accel_cam_orig.tFloat(tCropInds)-accel_cam_orig.tFloat(find(tCropInds,1)), 'measured',accel_UAV_dsamp.measured(tCropInds));
            gyro_UAV_dsamp_crop = struct('tFloat',gyro_UAV_dsamp.tFloat(tCropInds), 't',accel_cam_orig.tFloat(tCropInds)-accel_cam_orig.tFloat(find(tCropInds,1)), 'measured',gyro_UAV_dsamp.measured(tCropInds));

			% Save processed variables in output struct
            data(j).p_cam.(magnitudes{i}) = pos_cam_crop;
            data(j).v_cam.(magnitudes{i}) = vel_cam_crop;
            data(j).a_cam.(magnitudes{i}) = accel_cam_crop;
            data(j).a_UAV.(magnitudes{i}) = accel_UAV_dsamp_crop;
            data(j).gyro_UAV.(magnitudes{i}) = gyro_UAV_dsamp_crop;
			data(j).p_cam_orig.(magnitudes{i}) = pos_cam_orig;
			data(j).v_cam_orig.(magnitudes{i}) = vel_cam_orig;
            data(j).a_cam_orig.(magnitudes{i}) = accel_cam_orig;
            data(j).a_UAV_orig.(magnitudes{i}) = accel_UAV_crop;
            data(j).gyro_UAV_orig.(magnitudes{i}) = gyro_UAV_orig;
            data(j).tInterv = xl;
            data(j).tCropInds = tCropInds;
		end
		if ~isempty(drone_id_orig), data(j).drone_id = drone_id_crop; end
    end
end