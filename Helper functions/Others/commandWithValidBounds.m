%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensures a command has its components within a valid range: 0<=rho, 0<=theta<2*pi, 0<=phi<pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function commandOut = commandWithValidBounds(commandIn)
	commandIn = reshape(commandIn, 1,[]); % Make sure we treat it as a row vector (avoids errors later on)
	rho = reshape(commandIn(1:3:end), 1,[]); % Save rho, because whenever rho<0, we need to flip [multiply by -1] its corresponding theta and phi
	commandMatrix = [abs(rho); mod(sign(rho).*commandIn(2:3:end), 2*pi); pi/2 - (sign(rho).*(pi/2 - commandIn(3:3:end)))]; % 3xN matrix with rho, theta and phi in each row
	commandOut = reshape(commandMatrix, 1,[])'; % Reshape as [rho_1;theta_1;phi_1;rho_2;theta_2;phi_2;...]
end
