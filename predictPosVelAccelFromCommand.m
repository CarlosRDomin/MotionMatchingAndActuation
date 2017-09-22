%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicts what the acceleration, velocity and position of drone(s) would be given a polar-coord actuation command
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [accel, vel, pos] = predictPosVelAccelFromCommand(commandRho, commandTheta, commandPhi, spotterFps, deltaP, deltaT)
	persistent a v p t camFps dP dT;

	if isempty(a) || (nargin>=4 && camFps~=spotterFps) || (nargin>=5 && dP~=deltaP) || (nargin>=6 && dT~=deltaT)	% First time -> Compute a, v, p. Otherwise, no need to regenerate 1d synthetic motion
		if nargin<4 && isempty(camFps), spotterFps = 30; end
		if nargin<5 && isempty(dP), deltaP = 1; end
		if nargin<6 && isempty(dT), deltaT = 1; end
		camFps = spotterFps;
		dP = deltaP;
		dT = deltaT;
		
		K = 2*pi*deltaP/deltaT.^2; % a(t) = K*sin(2*pi*(1/deltaT)*t)
		t = 1/spotterFps : 1/spotterFps : deltaT;
		tt = 0 : 1/(10*spotterFps) : deltaT; % Avoid aliasing by generating the signal at a high sample rate

		aa = K*sin(2*pi/deltaT*tt); aa(end)=0; % Improve precision
		vv = cumtrapz(tt,aa); vv(end)=0; % Improve precision
		pp = cumtrapz(tt,vv);
		a = interp1(tt,aa, t);
		v = interp1(tt,vv, t);
		p = interp1(tt,pp, t);
	end

	[ax,ay,az] = sph2cart(commandTheta(:), pi/2-commandPhi(:), commandRho(:)*a);	% Use (:) so multiple signals can be computed at once by providing a vector of theta's, phi's and rho's (instead of a scalar)
	accel = [permute(ax,3:-1:1); permute(ay,3:-1:1); permute(az,3:-1:1)];	% To handle both scalar and vector theta/phi/rho inputs, we switch dimensions 1 and 3, so that ax,ay,az are always of size 1 x pointsPerIter x length(commandTheta) -> accel is 3 x pointsPerIter x length(commandTheta)
	accel(abs(accel)<10*eps)=0;	% Improve precision
	if nargout>1
		[vx,vy,vz] = sph2cart(commandTheta(:), pi/2-commandPhi(:), commandRho(:).*v); vel = [permute(vx,3:-1:1);permute(vy,3:-1:1);permute(vz,3:-1:1)]; vel(abs(vel)<10*eps)=0;
	end
	if nargout>2
		[px,py,pz] = sph2cart(commandTheta(:), pi/2-commandPhi(:), commandRho(:).*p); pos = [permute(px,3:-1:1);permute(py,3:-1:1);permute(pz,3:-1:1)]; pos(abs(pos)<10*eps)=0;
	end
end
