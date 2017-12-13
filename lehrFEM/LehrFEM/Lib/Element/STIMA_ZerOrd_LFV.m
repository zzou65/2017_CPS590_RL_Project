function aLoc = STIMA_ZerOrd_LFE(vertices,midPoints,center,method,bdFlags,rHandle,varargin)
% STIMA_ZERORD_LFV Element stiffness matrix for the zero-order term.
%
%   ALOC = STIMA_GENLAPL_LFV(VERTICES,MIDPOINTS,CENTER,KHANDLE) computes the 
%   element stiffness matrix for the zero-order term  using linear
%   finite volumes.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   MIDPOINTS is 3-by-2 matrix specifying the border points on the edges of
%   the element.
%
%   CENTER is a length 2 row vector specifying the center point of the element.
%
%   CHANDLE is a function handle for the r-function.
%
%   Example:
%
%   Aloc = STIMA_ZerOrd_LFV([0 0; 1 0; 0 1],@(x)1);
%
%   Copyright 2007-2007 Eivind Fonn
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

	% Compute required coefficients

	m(1) = (abs(det([vertices(1,:) 1;midPoints(1,:) 1;center 1])) + abs(det([vertices(1,:) 1;midPoints(3,:) 1;center 1])))/2;
	m(2) = (abs(det([vertices(2,:) 1;midPoints(1,:) 1;center 1])) + abs(det([vertices(2,:) 1;midPoints(2,:) 1;center 1])))/2;
	m(3) = (abs(det([vertices(3,:) 1;midPoints(2,:) 1;center 1])) + abs(det([vertices(3,:) 1;midPoints(3,:) 1;center 1])))/2;

	r(1) = rHandle(vertices(1,:));
	r(2) = rHandle(vertices(2,:));
	r(3) = rHandle(vertices(3,:));

	% Compute stiffness matrix

	aLoc = zeros(3,3);

	for i=1:3
		aLoc(i,i) = r(i)*m(i);
	end
return
