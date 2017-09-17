%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a test. Tries to communicate with the python CF script and send it the loop info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateScoreLoop
    s = tcpip('0.0.0.0', 30000, 'NetworkRole', 'server', 'Timeout', 30.0);
    s.BytesAvailableFcnMode = 'terminator';
    s.Terminator = 10;
    s.BytesAvailableFcn = @onData;
    fopen(s);
    inited = false;

    N=0; M=0;
    % get the N and M values

    ready = false;
    dataLen = 0;
    nRec = 0;
    rPrior = 0;
    pos = 0;
    function onData(sock, ~)
        if ~inited
            allData = fscanf(sock, '%f');
            N = allData(1); M = allData(2);
            dataLen = N*M+N*3;
            inited = true;
        else
            allData = fscanf(sock, '%f');
            rPrior = reshape(allData(1:N*M), N, M);
            pos = reshape(allData(N*M+1:dataLen), N,3);
            nRec = nRec + 1;
            ready = true;
        end
    end
    while nRec < 5
        pause(0.5);
        if ready
           comm = doLoopCommand(rPrior, pos, N);
           fprintf(s, '%f,', comm);
           ready = false;
        end
    end
    fclose(s);
end

function [commandOut] = doLoopCommand(runningPrior, posUAVcam, N)
    P_total = 0;
    expectedCommand = zeros(3*N,1);
    deltaP = 1;	% 1m per timestep
    deltaT = 1;	% Timestep/iteration: 1s
    minRisk = 0.5; % m, min closest distance allowed
    sigmaNoiseCam = 0.05;	% m, noise in camera's position estimation
    sigmaNoisePos = 0.015;	% m, error in position after performing a motion command
    sigmaNoiseAccel = 0.25;	% m/s2, noise in IMU's accelerometer data
    
    nAssignments = 3;
    assignmentList = computeNBestAssignments(nAssignments, runningPrior, -log(1e-30), -log(1e-2));

    for i = 1:nAssignments
        assignment = assignmentList(i).matches;
        P_assignment = prod(runningPrior(sub2ind(size(runningPrior), assignment(:,1),assignment(:,2))));
        P_total = P_total + P_assignment;
        %sortedAssignment = sortrows([assignments; [unassignedUAVs(:) NaN(length(unassignedUAVs),1)]]);
        [command,newP,exitFlag,output] = fmincon(@(x) (-estimateImprovementOfCommand(x,assignment,runningPrior,[],sigmaNoiseAccel)), ...
            zeros(3*N,1),[],[],[],[],zeros(3*N,1),repmat([deltaP,2*pi,pi]',N,1), ...
            @(x) deal(minRisk-estimateRiskOfCommand(x,assignment,posUAVcam), 0)); %, optimoptions('fmincon', 'Algorithm','active-set'));
        expectedCommand = expectedCommand + P_assignment.*command;
    end
    commandOut = expectedCommand./P_total;
    dispImproved(sprintf('Sending motion command:\nrho:\t%s\ntheta:\t%s\nphi:\t%s\n', num2str(commandOut(1:3:end)','%8.2f'), num2str(rad2deg(commandOut(2:3:end))','%7.1f?'), num2str(rad2deg(commandOut(3:3:end))','%7.1f?')), 'keepthis');

end

