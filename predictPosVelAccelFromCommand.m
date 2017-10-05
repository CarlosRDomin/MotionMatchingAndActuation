%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicts what the acceleration, velocity and position of drone(s) would be given a polar-coord actuation command
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [accel, vel, pos] = predictPosVelAccelFromCommand(commandRho, commandTheta, commandPhi, a0, v0, p0, spotterFps, deltaT, sigmaEnvironmentNoise)
	if nargin<4 || isempty(a0), a0 = 0; end
	if nargin<5 || isempty(v0), v0 = 0; end
	if nargin<6 || isempty(p0), p0 = 0; end
	if nargin<7 || isempty(spotterFps), spotterFps = 30; end
	if nargin<8 || isempty(deltaT), deltaT = 1; end
	if nargin<9 || isempty(sigmaEnvironmentNoise), sigmaEnvironmentNoise = 0.75; end
	
	% Compute how much we want to move in each cartesian direction
	commandRho = reshape(commandRho, [],1);
	commandTheta = reshape(commandTheta, [],1);
	commandElev = reshape(pi/2 - commandPhi, [],1);
	[dx,dy,dz] = sph2cart(commandTheta, commandElev, commandRho); % Convert spherical to cartesian coords
	deltaP = cat(3, dx, dy, dz); % numel(commandRho) x 1 x 3

	t = 1/spotterFps : 1/spotterFps : deltaT;
	tt = 0 : 1/(10*spotterFps) : deltaT; % Avoid aliasing by generating the signal at a high sample rate
	
	% Compute constants for our command model: a(t) = K*sin(2*pi*(1/deltaT)*t) + At^3 + Bt^2 + Ct + a0 [+ noiseA]
		% syms a0 v0 p0 deltaT deltaP A B C K
		% K = 2*pi.*deltaP./deltaT.^2;
		% eq1 = deltaP == deltaT.^2/2*(K./pi+a0) + deltaT.*v0 + deltaT.^5/20.*A + deltaT.^4/12.*B + deltaT.^3/6.*C;
		% eq2 = 0 == deltaT.*a0 + v0 + deltaT.^4/4.*A + deltaT.^3/3.*B + deltaT.^2/2.*C;
		% eq3 = 0 == a0 + deltaT.^3.*A + deltaT.^2.*B + deltaT.*C;
		% sol = solve([eq1 eq2 eq3], [A B C]);
		% A = sol.A; B = sol.B; C = sol.C;
	K = 2*pi.*deltaP./deltaT.^2;
	A = -(10*(6*deltaT^2.*K - 12*pi.*deltaP + 6*pi*deltaT.*v0 + pi*deltaT^2.*a0))./(deltaT^5*pi);
	B = (6*(15*deltaT^2.*K - 30*pi.*deltaP + 16*pi*deltaT.*v0 + 3*pi*deltaT^2.*a0))./(deltaT^4*pi);
	C = -(3*(10*deltaT^2.*K - 20*pi.*deltaP + 12*pi*deltaT.*v0 + 3*pi*deltaT^2.*a0))./(deltaT^3*pi);
	noiseA = cat(2, zeros(numel(commandRho),1,3), cumsum(sigmaEnvironmentNoise/sqrt(numel(tt)-1) .* randn(numel(commandRho), numel(tt)-1, 3), 2)); % Random walk on aa. Sigma of the distance after n points is sqrt(n)*sigmaOrig -> ~68% of the points will be within 0±sqrt(n)*sigmaOrig -> sigmaEnvironmentNoise = sqrt(numel(tt))*sigmaOrig -> sigmaOrig = sigmaEnvironmentNoise/sqrt(numel(tt))
	
	% Compute acceleration and decimate to the actual sampling rate
	aa = a0 + K.*sin(2*pi/deltaT.*tt) + A.*tt.^3 + B.*tt.^2 + C.*tt + noiseA;
	accel = permute(interp1(tt,permute(aa, [2 1 3]), t), [2 1 3]);
	
	% If needed, integrate to compute vel and pos
	if nargout > 1
		vv = v0 + cumtrapz(tt, aa, 2);
		vel = permute(interp1(tt,permute(vv, [2 1 3]), t), [2 1 3]);
	end
	if nargout > 2
		pp = p0 + cumtrapz(tt, vv, 2);
		pos = permute(interp1(tt,permute(pp, [2 1 3]), t), [2 1 3]);
	end
end
