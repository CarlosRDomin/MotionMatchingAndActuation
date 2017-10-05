function signalOut = derivFilter(signalIn, deriv, dt, poly_order, win_size)
	% signalIn is a Nxlength(t)xdims, where N indicates independent signals (eg from different drones), length(t) is the time duration, dims is the signal dimensions (eg 3D)
	if nargin<4 || isempty(poly_order)
		poly_order = 2;
	end
	if nargin<5 || isempty(win_size)
		win_size = 21;
	end
	win_half_size = (win_size-1)/2;
	
% 	[~, diffMatrix] = sgolay(poly_order, win_size);
% 	diffVect = factorial(deriv)/(-dt)^deriv .* diffMatrix(:,deriv+1);
	diffMatrix = sgDerivFilterCoefs(dt*(-win_half_size:win_half_size), poly_order,deriv);
	signalOut = zeros(size(signalIn));
	for iSignal = 1:size(signalIn,1)
		for iDim = 1:size(signalIn,3)
			% Compute steady state
			signalOut(iSignal,:,iDim) = conv(signalIn(iSignal,:,iDim), diffMatrix(win_half_size+1,:), 'same');
			
			% Compute initial transient
			beginningIndsIn = 1 : min(win_size, size(signalIn,2));
			beginningIndsOut = 1 : min(win_half_size, size(signalIn,2));
			signalOut(iSignal,beginningIndsOut,iDim) = reshape(diffMatrix(beginningIndsOut,beginningIndsIn)*reshape(signalIn(iSignal,beginningIndsIn,iDim), [],1), 1,[]);
			
			% Compute ending transient
			endingIndsIn = size(signalIn,2)+1-beginningIndsIn(end:-1:1);
			endingIndsOut = size(signalIn,2)+1-beginningIndsOut(end:-1:1);
			signalOut(iSignal,endingIndsOut,iDim) = reshape(diffMatrix(end-numel(endingIndsOut)+1:end,end-numel(endingIndsIn)+1:end)*reshape(signalIn(iSignal,endingIndsIn,iDim), [],1), 1,[]);
		end
	end
end
