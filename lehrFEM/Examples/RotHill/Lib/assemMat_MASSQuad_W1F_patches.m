function varargout= assemMat_MASSQuad_W1F_patches(Mesh, defMesh,intersec,varargin)
% 
% Mesh structure; tracedvertices array with coordinates of traced vertices and
% local local elements; Jac function handle to Jacobian of velocity field,
% h timestepsize
%
%   Copyright 2008-2009 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

% compute the total number of patches
nPatches = 0;
for i = 1:nElements
    nPatches = nPatches + intersec(i).nElems;
end

% Preallocate memory

I = zeros(9*nPatches,1);
J = zeros(9*nPatches,1);
A = zeros(9*nPatches,1);

% Check for element flags
if (isfield(Mesh,'ElemFlag')), flags = Mesh.ElemFlag;
else flags = zeros(nElements,1); end

% Assemble element contributions
loc = 1:9;
for i = 1:nElements
  polygons = patches(Mesh,defMesh,i,intersec(i));
  for j = 1:intersec(i).nElems  
    Aloc=zeros(3,3);
    
    % Vertices
    vid = Mesh.Elements(i,:);
    a1 = Mesh.Coordinates(vid(1),:);
    a2 = Mesh.Coordinates(vid(2),:);
    a3 = Mesh.Coordinates(vid(3),:);
 
    % Compute element mapping
    bK = a1;
    BK = [a2-bK; ...
        a3-bK];
    det_BK = abs(det(BK));
    inv_BK = inv(BK);
    TK = transpose(inv_BK);
    
    % compute barycenter of the intersection
    p_Element = polygons(j).Elem;
    polygon = polygons(j).Polygon;
    evalpoint = sum(polygon,1)/size(polygon,1);
   %evalpoint = sum(polygon(1:2,:),1)/2;
    
    % transfer the evalutation point to the standard triangle
    stevalpoint = (evalpoint - a1)*inv_BK;

    % Extract global edge numbers

    eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
        Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
        Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];
 
    N = shap_W1F(stevalpoint);
    N(:,[1 2]) = N(:,[1 2])*TK;
    N(:,[3 4]) = N(:,[3 4])*TK;
    N(:,[5 6]) = N(:,[5 6])*TK;
    
     % Determine the orientation
    if(Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else    p1 = -1;  end 
    if(Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else    p3 = -1;  end
    
    % Add contributions to stiffness matrix
    Aloc = ...
        [p1*N([1 2])*N([1 2])'*p1 p1*N([1 2])*N([3 4])'*p2 p1*N([1 2])*N([5 6])'*p3; ...
        p2*N([3 4])*N([1 2])'*p1 p2*N([3 4])*N([3 4])'*p2 p2*N([3 4])*N([5 6])'*p3; ...
        p3*N([5 6])*N([1 2])'*p1 p3*N([5 6])*N([3 4])'*p2 p3*N([5 6])*N([5 6])'*p3];
        
    Aloc = Aloc*polyarea(polygon(:,1),polygon(:,2));
    
    I(loc) = set_Rows(eidx,3);
    J(loc) = set_Cols(eidx,3);
    A(loc) = Aloc(:);
    loc = loc+9;
  end 
end

% Assign output arguments
if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I,J,A,nEdges,nEdges);
end