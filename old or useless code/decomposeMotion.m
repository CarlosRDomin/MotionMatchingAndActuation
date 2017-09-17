function motionMagns = decomposeMotion(x, n, base)
	if nargin<3 || isempty(base)
		base = 3;
	end

	motionMagns = repmat(-1, 1,n);
	aux = dec2base(x-1, base);
	for i = 1:length(aux)
		motionMagns(i) = str2double(aux(length(aux)-i+1))-1;
% 		motionMagns(i) = mod(x, base);
% 		x = floor(x./base);
	end
end
