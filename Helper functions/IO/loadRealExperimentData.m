%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads (real, ie. from CF drones) experiment data from .npz files into a struct which contains the accel measured from the camera, accel and gyro measured from the drone's IMU...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = loadRealExperimentData(logInfo, logFolder)
    magnitudes = {'X'; 'Y'; 'Z'};
    if nargin < 2 || isempty(logFolder)
        logFolder = [fileparts(mfilename('fullpath')) '/../../data/Real'];
    end
    
    data = repmat(cell2struct([repmat({cell2struct(cell(3,1), magnitudes)}, 6,1); cell(2,1)], {'a_cam','a_UAV','gyro_UAV','a_cam_orig','a_UAV_orig','gyro_UAV_orig', 'tInterv', 'tCropInds'}), length({logInfo.datetime}), 1);
    for j = 1:length({logInfo.datetime})
        for i = 1:size(magnitudes,1)
            accel_cam_orig = readNPZ([logFolder '/' logInfo(j).datetime '/log_a' magnitudes{i} '_cam_' logInfo(j).ch '_' strrep(logInfo(j).datetime, ' ','_') '.npz']);
            accel_UAV_orig = readNPZ([logFolder '/' logInfo(j).datetime '/log_a' magnitudes{i} '_world_' logInfo(j).ch '_' strrep(logInfo(j).datetime, ' ','_') '.npz']);
            gyro_UAV_orig = readNPZ([logFolder '/' logInfo(j).datetime '/log_gyro' magnitudes{i} '_' logInfo(j).ch '_' strrep(logInfo(j).datetime, ' ','_') '.npz']);
            accel_UAV_dsamp = struct('tFloat',accel_cam_orig.tFloat, 'measured',interp1(accel_UAV_orig.tFloat, 9.81*movingAvgFilter(30, accel_UAV_orig.measured), accel_cam_orig.tFloat)); % accel_UAV_dsamp.measured(isnan(accel_UAV_dsamp.measured))=0;
            gyro_UAV_dsamp = struct('tFloat',accel_cam_orig.tFloat, 'measured',interp1(gyro_UAV_orig.tFloat, movingAvgFilter(30, gyro_UAV_orig.measured), accel_cam_orig.tFloat)); % gyro_UAV_dsamp.measured(isnan(gyro_UAV_dsamp.measured))=0;

            xl = [accel_cam_orig.tFloat(6), accel_UAV_orig.tFloat(end)];
            tCropInds = (accel_cam_orig.tFloat>=xl(1) & accel_cam_orig.tFloat<=xl(2));
            accel_cam_crop = struct('tFloat',accel_cam_orig.tFloat(tCropInds), 't',accel_cam_orig.tFloat(tCropInds)-accel_cam_orig.tFloat(find(tCropInds,1)) , 'measured',accel_cam_orig.measured(tCropInds));
            accel_UAV_dsamp_crop = struct('tFloat',accel_UAV_dsamp.tFloat(tCropInds), 't',accel_cam_orig.tFloat(tCropInds)-accel_cam_orig.tFloat(find(tCropInds,1)), 'measured',accel_UAV_dsamp.measured(tCropInds));
            gyro_UAV_dsamp_crop = struct('tFloat',gyro_UAV_dsamp.tFloat(tCropInds), 't',accel_cam_orig.tFloat(tCropInds)-accel_cam_orig.tFloat(find(tCropInds,1)), 'measured',gyro_UAV_dsamp.measured(tCropInds));

            data(j).a_cam.(magnitudes{i}) = accel_cam_crop;
            data(j).a_UAV.(magnitudes{i}) = accel_UAV_dsamp_crop;
            data(j).gyro_UAV.(magnitudes{i}) = gyro_UAV_dsamp_crop;
            data(j).a_cam_orig.(magnitudes{i}) = accel_cam_orig;
            data(j).a_UAV_orig.(magnitudes{i}) = accel_UAV_orig;
            data(j).gyro_UAV_orig.(magnitudes{i}) = gyro_UAV_orig;
            data(j).tInterv = xl;
            data(j).tCropInds = tCropInds;
        end
    end
end