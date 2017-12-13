function varargout= assemMat_SemiLagQuad_W1F_bary(Mesh,tracedbcenter,varargin)
% assemMat_SemiLagQuad_W1F_bary calculates the stiffness matrix for 
% Semi-Lagrange for Whitney one forms using simple vertex based
% quadrature rule
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

% Preallocate memory
nElements = size(Mesh.Elements,1);

% Preallocate memory

I = zeros(9*nElements,1);
J = zeros(9*nElements,1);
A = zeros(9*nElements,1);

% Check for element flags
if (isfield(Mesh,'ElemFlag')), flags = Mesh.ElemFlag;
else flags = zeros(nElements,1); end

% Assemble element contributions
loc = 1:9;
for i = 1:nElements
    
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

    % Extract global edge numbers

    eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
        Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
        Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];
 
    N = shap_W1F([1/3 1/3]);
    N(:,[1 2]) = N(:,[1 2])*TK;
    N(:,[3 4]) = N(:,[3 4])*TK;
    N(:,[5 6]) = N(:,[5 6])*TK;
    
    point=tracedbcenter(i,[1,2]);
    p_Element=tracedbcenter(i,3);

    % Vertices
    p_vid = Mesh.Elements(p_Element,:);
    p_a1 = Mesh.Coordinates(p_vid(1),:);
    p_a2 = Mesh.Coordinates(p_vid(2),:);
    p_a3 = Mesh.Coordinates(p_vid(3),:);
    
   p_eidx= [Mesh.Vert2Edge(Mesh.Elements(p_Element,2),Mesh.Elements(p_Element,3)) ...
        Mesh.Vert2Edge(Mesh.Elements(p_Element,3),Mesh.Elements(p_Element,1)) ...
        Mesh.Vert2Edge(Mesh.Elements(p_Element,1),Mesh.Elements(p_Element,2))];
    
    % Determine the orientation

    if(Mesh.Edges(p_eidx(1),1)==p_vid(2)),  p_p1 = 1;  else    p_p1 = -1;  end
    if(Mesh.Edges(p_eidx(2),1)==p_vid(3)),  p_p2 = 1;  else    p_p2 = -1;  end
    if(Mesh.Edges(p_eidx(3),1)==p_vid(1)),  p_p3 = 1;  else    p_p3 = -1;  end
    
    % Compute element mapping
    p_bK= p_a1;
    p_BK = [p_a2-p_bK; ...
        p_a3-p_bK];
    p_det_BK = abs(det(p_BK));
    p_inv_BK = inv(p_BK);
    p_TK = transpose(p_inv_BK);
    
    point_hat=(point-p_bK)*p_inv_BK;
    p_N =shap_W1F(point_hat);
    
    WM=[tracedbcenter(i,4) tracedbcenter(i,6);...
        tracedbcenter(i,5) tracedbcenter(i,7)];
    
    p_N([1 2]) = p_N([1 2])*p_TK*WM;
    p_N([3 4]) = p_N([3 4])*p_TK*WM;
    p_N([5 6]) = p_N([5 6])*p_TK*WM;
    
    % Determine the orientation
    if(Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else    p1 = -1;  end
    if(Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else    p3 = -1;  end
    
    % Add contributions to stiffness matrix
    Aloc = ...
        [p1*N([1 2])*p_N([1 2])'*p_p1 p1*N([1 2])*p_N([3 4])'*p_p2 p1*N([1 2])*p_N([5 6])'*p_p3; ...
        p2*N([3 4])*p_N([1 2])'*p_p1 p2*N([3 4])*p_N([3 4])'*p_p2 p2*N([3 4])*p_N([5 6])'*p_p3; ...
        p3*N([5 6])*p_N([1 2])'*p_p1 p3*N([5 6])*p_N([3 4])'*p_p2 p3*N([5 6])*p_N([5 6])'*p_p3];
        
    Aloc = Aloc*det_BK/2;
    
    I(loc) = set_Rows(eidx,3);    
    J(loc) = set_Cols(p_eidx,3);
    A(loc) = Aloc(:);
    loc = loc+9;
    
end

% Assign output arguments
if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I,J,A,nEdges,nEdges);
end
