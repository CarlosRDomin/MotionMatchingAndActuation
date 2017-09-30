function [points, lambda] = intersectLineSphere(line, sphere, tolerance)
%INTERSECTLINESPHERE Return intersection points between a line and a sphere
% https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d?focused=7621190&tab=function
%
%   PTS = intersectLineSphere(LINE, SPHERE);
%   Returns the two points which are the intersection of the given line and
%   sphere. 
%   LINE   : [x0 y0 z0  dx dy dz]
%   SPHERE : [xc yc zc  R]
%   PTS     : [x1 y1 z1 ; x2 y2 z2]
%   If there is no intersection between the line and the sphere, return a
%   2-by-3 array containing only NaN.
%
%   Example
%     % draw the intersection between a sphere and a collection of parallel
%     % lines 
%     sphere = [50.12 50.23 50.34 40];
%     [x, y] = meshgrid(10:10:90, 10:10:90);
%     n = numel(x);
%     lines = [x(:) y(:) zeros(n,1) zeros(n,2) ones(n,1)];
%     figure; hold on; axis equal;
%     axis([0 100 0 100 0 100]); view(3);
%     drawSphere(sphere);
%     drawLine3d(lines);
%     pts = intersectLineSphere(lines, sphere);
%     drawPoint3d(pts, 'rx');
%
%     % apply rotation on set of lines to check with non vertical lines
%     rot = eulerAnglesToRotation3d(20, 30, 10);
%     rot2 = recenterTransform3d(rot, [50 50 50]);
%     lines2 = transformLine3d(lines, rot2);
%     figure; hold on; axis equal;
%     axis([0 100 0 100 0 100]); view(3);
%     drawSphere(sphere);
%     drawLine3d(lines2);
%     pts2 = intersectLineSphere(lines2, sphere);
%     drawPoint3d(pts, 'rx');
%
%   See also
%   spheres, circles3d, intersectPlaneSphere
%

%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 18/02/2005.
%

%   HISTORY
%   2011-06-21 bug for tangent lines, add tolerance


	%% Process input arguments
	
	% Extend line or sphere vertically if needed (test 1-to-many)
	if size(line,1)~=size(sphere,1)
		if min(size(line,1), size(sphere,1)) > 1 % N lines and M spheres not allowed 
			error('MatGeom:geom3d:intersectLineSphere', 'Inputs must have the same number of rows, or one must have size 1 in their first dimension');
		else
			if size(line,1)<size(sphere,1), line=repmat(line, size(sphere,1),1); else, sphere = repmat(sphere, size(line,1),1); end
		end
	end
	if nargin<3 || isempty(tolerance), tolerance = 1e-14; end

	% difference between centers
	dc = bsxfun(@minus, line(:, 1:3), sphere(:, 1:3));

	% equation coefficients
	a = sum(line(:, 4:6) .* line(:, 4:6), 2);
	b = 2 * sum(bsxfun(@times, dc, line(:, 4:6)), 2);
	c = sum(dc.*dc, 2) - sphere(:,4).*sphere(:,4);

	% solve equation
	delta = b.*b - 4*a.*c;

	% initialize empty results
	points = NaN(2*size(delta,1), 3);
	lambda = NaN(2*size(delta,1), 1);


	%% Process intersection points

	% process couples with two intersection points
	delta(abs(delta)<tolerance) = 0;
	inds = find(delta >= 0);
	if ~isempty(inds)
		% delta positive: find two roots of second order equation
		lambda(2*inds-1) = (-b(inds) -sqrt(delta(inds))) / 2 ./ a(inds);
		lambda(2*inds) = (-b(inds) +sqrt(delta(inds))) / 2 ./ a(inds);

		% convert into 3D coordinate
		points(2*inds-1, :) = line(inds, 1:3) + bsxfun(@times, lambda(2*inds-1), line(inds, 4:6));
		points(2*inds,:) = line(inds, 1:3) + bsxfun(@times, lambda(2*inds), line(inds, 4:6));
	end
end