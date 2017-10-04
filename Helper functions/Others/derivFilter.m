function signalOut = derivFilter(signalIn, deriv, dt, poly_order, win_size)
	% signalIn is a Nxlength(t)xdims, where N indicates independent signals (eg from different drones), length(t) is the time duration, dims is the signal dimensions (eg 3D)
	if nargin<4 || isempty(poly_order)
		poly_order = 2;
	end
	if nargin<5 || isempty(win_size)
		win_size = 21;
	end
	%win_half_size = (win_size-1)/2;
	
	[~, diffMatrix] = sgolay(poly_order, win_size);
	diffVect = factorial(deriv)/(-dt)^deriv .* diffMatrix(:,deriv+1);
	signalOut = zeros(size(signalIn));
	for iSignal = 1:size(signalIn,1)
		for iDim = 1:size(signalIn,3)
			signalOut(iSignal,:,iDim) = conv(signalIn(iSignal,:,iDim), diffVect, 'same');
			%signalOut(iSignal,1:win_half_size,iDim) = reshape(repmat(diffVect', win_half_size,1)*reshape(signalIn(iSignal,1:win_size,iDim), [],1), 1,[]);
			%signalOut(iSignal,1:win_half_size,iDim) = reshape(repmat(diffVect', win_half_size,1)*reshape(signalIn(iSignal,end-win_size+1:end,iDim), [],1), 1,[]);
		end
	end
end
