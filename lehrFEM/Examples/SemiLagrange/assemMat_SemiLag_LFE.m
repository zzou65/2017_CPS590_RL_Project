function varargout= assemMat_SemiLag_LFE(Mesh, tracedvertices, varargin)
% assemMat_SemiLag_LFE calculates the interpolation matrix for 
% Semi-Lagrange for linear finite elements.
%
% Mesh structure; tracedvertices array with coodinates of traced vertices and
% local local elements.
%
% Example:
% M_0 =  assemMat_Mass0fD(NewMesh);
% pulled_back_vertices=trace_vertices(NewMesh,-dt*directions)
% A =assemMat_SemiLag_LFE(Mesh,pulled_back_vertices);
% B =assemMat_SemiLag_Quad_LFE(Mesh,pulled_back_vertices); 
% % M_0 * A == B;
%
% Copyright 2008-2009 Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

% Preallocate memory
A = zeros(3*nCoordinates,1);
I = zeros(3*nCoordinates,1);
J = zeros(3*nCoordinates,1);

loc = [1 2 3];

% loop over all elements and (back)-trace vertices
for i = 1:nCoordinates
    
    point = tracedvertices(i,[1,2]);
    Element = tracedvertices(i,3);
   
    % Vertices
    vid = Mesh.Elements(Element,:);
    a1 = Mesh.Coordinates(vid(1),:);
    a2 = Mesh.Coordinates(vid(2),:);
    a3 = Mesh.Coordinates(vid(3),:);
    
    % Compute element mapping
    bK = a1;
    BK = [a2-bK; ...
        a3-bK];
    det_BK = abs(det(BK));
    inv_BK = inv(BK);
    
    point_hat=(point-bK)*inv_BK;
    
    A(loc(1)) = 1-sum(point_hat);
    A(loc(2)) = point_hat(1);
    A(loc(3)) = point_hat(2);
    I(loc) = [i i i];
    J(loc) = vid;
    loc=loc+3;
end

if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I,J,A,nCoordinates,nCoordinates);;
end
