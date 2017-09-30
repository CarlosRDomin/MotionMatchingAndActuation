function projVect = projectVectPlane(vect, planeNormalVect)
	% proj_plane(u) = u - proj_n(u) = u - (u•n)/(n•n)*n
	projVect = vect - dot(vect, planeNormalVect)/dot(planeNormalVect, planeNormalVect) * planeNormalVect;
end
