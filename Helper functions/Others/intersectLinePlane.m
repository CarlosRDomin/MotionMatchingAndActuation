function [point, lambda] = intersectLinePlane(line, plane, tolerance)
%INTERSECTLINEPLANE Intersection point between a 3D line and a plane
% https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d?focused=7621187&tab=function
%
%   PT = intersectLinePlane(LINE, PLANE)
%   Returns the intersection point of the given line and the given plane.
%   LINE:  [x0 y0 z0 dx dy dz]
%   PLANE: [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
%   PT:    [xi yi zi]
%   If LINE and PLANE are parallel, return [NaN NaN NaN].
%   If LINE (or PLANE) is a matrix with 6 (or 9) columns and N rows, result
%   is an array of points with N rows and 3 columns.
%   
%   PT = intersectLinePlane(LINE, PLANE, TOL)
%   Specifies the tolerance factor to test if a line is parallel to a
%   plane. Default is 1e-14.
%
%   Example
%     % define horizontal plane through origin
%     plane = [0 0 0   1 0 0   0 1 0];
%     % intersection with a vertical line
%     line = [2 3 4  0 0 1];
%     intersectLinePlane(line, plane)
%     ans = 
%        2   3   0
%     % intersection with a line "parallel" to plane
%     line = [2 3 4  1 2 0];
%     intersectLinePlane(line, plane)
%     ans = 
%       NaN  NaN  NaN
%
%   See also:
%   lines3d, planes3d, points3d, clipLine3d
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 17/02/2005.
%

%   HISTORY
%   24/11/2005 add support for multiple input
%   23/06/2006 correction from Songbai Ji allowing different number of
%       lines or plane if other input has one row
%   14/12/2006 correction for parallel lines and plane normals
%   05/01/2007 fixup for parallel lines and plane normals
%   24/04/2007 rename as 'intersectLinePlane'
%   11/19/2010 Added bsxfun functionality for improved speed (Sven Holcombe)
%   01/02/2011 code cleanup, add option for tolerance, update doc


	%% Process input arguments
	
	% Extend line or sphere vertically if needed (test 1-to-many)
	if size(line,1)~=size(plane,1) && min(size(line,1), size(plane,1)) > 1 % N planes and M lines not allowed 
		error('MatGeom:geom3d:intersectLinePlane', 'Inputs must have the same number of rows, or one must have size 1 in their first dimension');
	end
	if nargin<3 || isempty(tolerance), tolerance = 1e-14; end

	% plane normal
	n = cross(plane(:,4:6), plane(:,7:9), 2);

	% difference between origins of plane and line
	dp = bsxfun(@minus, plane(:, 1:3), line(:, 1:3));

	% dot product of line direction with plane normal
	denom = sum(bsxfun(@times, n, line(:,4:6)), 2);

	% relative position of intersection point on line (can be inf in case of a
	% line parallel to the plane)
	lambda = sum(bsxfun(@times, n, dp),2) ./ denom;
	lambda(abs(denom)<tolerance) = NaN;

	% compute coord of intersection point
	point = bsxfun(@plus, line(:,1:3),  bsxfun(@times, repmat(lambda, 1,3), line(:,4:6)));
end