function L = assemLoad_LFV(Mesh,fHandle,varargin)
% ASSEMLOAD_LFV Assemble linear FV contributions.
%
%   L = ASSEMLOAD_LFV(MESH,FHANDLE) assembles the global load 
%   vector for the load data given by the function handle FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%    EDGES        L-by-2 matrix specifying the edges of the mesh.
%    MIDPOINTS    L-by-2 matrix specifying the mid-edge borders of the mesh.
%    CENTERPOINTS N-by-2 matrix specifying the mid-element points of the mesh.
%
%   L = ASSEMLOAD_LFE(MESH,FHANDLE,FPARAM) also handles the 
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)x(:,1).^2+x(:,2).^2;
%   L = assemLoad_LFV(Mesh,FHandle);

%   Copyright 2007-2007 Eivind Fonn
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constants
  
    nCoordinates = size(Mesh.Coordinates,1);
    nElements = size(Mesh.Elements,1);
  
    % Preallocate memory
  
    L = zeros(nCoordinates,1);
  
    % Compute load vector
  
    for i = 1:nElements
	vxid = Mesh.Elements(i,:);
	vertices = Mesh.Coordinates(vxid,:);
	edges = [Mesh.Vert2Edge(vxid(1),vxid(2));...
		 Mesh.Vert2Edge(vxid(2),vxid(3));...
		 Mesh.Vert2Edge(vxid(3),vxid(1))];
	midpoints = Mesh.MidPoints(edges,:);
	center = Mesh.CenterPoints(i,:);

	m(1) = (abs(det([vertices(1,:) 1;midpoints(1,:) 1;center 1])) + abs(det([vertices(1,:) 1;midpoints(3,:) 1;center 1])))/2;
	m(2) = (abs(det([vertices(2,:) 1;midpoints(1,:) 1;center 1])) + abs(det([vertices(2,:) 1;midpoints(2,:) 1;center 1])))/2;
	m(3) = (abs(det([vertices(3,:) 1;midpoints(2,:) 1;center 1])) + abs(det([vertices(3,:) 1;midpoints(3,:) 1;center 1])))/2;

	L(vxid(1)) = L(vxid(1)) + m(1)*fHandle(vertices(1,:),varargin);
	L(vxid(2)) = L(vxid(2)) + m(2)*fHandle(vertices(2,:),varargin);
	L(vxid(3)) = L(vxid(3)) + m(3)*fHandle(vertices(3,:),varargin);
    end
return
