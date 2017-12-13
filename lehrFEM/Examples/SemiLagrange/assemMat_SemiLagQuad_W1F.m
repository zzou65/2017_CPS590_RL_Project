function varargout= assemMat_SemiLagQuad_W1F(Mesh,tracedvertices,varargin)
% assemMat_SemiLagQuad_W1F calculates the stiffness matrix for 
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
I = zeros(27*nElements,1);
J = zeros(27*nElements,1);
A = zeros(27*nElements,1);

% Check for element flags
if (isfield(Mesh,'ElemFlag')), flags = Mesh.ElemFlag;
else flags = zeros(nElements,1); end

% Assemble element contributions
loc = 1:9;
for i = 1:nElements
    
    Aloc1=zeros(3,3);
    Aloc2=zeros(3,3);
    Aloc3=zeros(3,3);

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

    % loop over all three vertices
    Map_vid_eidx = zeros(3,3);
    Map_vid_values = zeros(3,6);
    Map_vid_orient = zeros(3,3);
    
    N = shap_W1F([0 0; 1 0; 0 1]);
    N(:,[1 2]) = N(:,[1 2])*TK;
    N(:,[3 4]) = N(:,[3 4])*TK;
    N(:,[5 6]) = N(:,[5 6])*TK;
        
    for j=1:3
        point = tracedvertices(vid(j),[1,2]);
        p_Element = tracedvertices(vid(j),3);
        
        WM = [tracedvertices(vid(j),4) tracedvertices(vid(j),6) ;...
        tracedvertices(vid(j),5) tracedvertices(vid(j),7)];
        
        % Vertices
        p_vid = Mesh.Elements(p_Element,:);
        p_a1 = Mesh.Coordinates(p_vid(1),:);
        p_a2 = Mesh.Coordinates(p_vid(2),:);
        p_a3 = Mesh.Coordinates(p_vid(3),:);
        
        Map_vid_eidx(j,:)= ...
            [Mesh.Vert2Edge(p_vid(2), p_vid(3)) ...
             Mesh.Vert2Edge(p_vid(3), p_vid(1)) ...
             Mesh.Vert2Edge(p_vid(1), p_vid(2))];

        % Determine the orientation
        if(Mesh.Edges(Map_vid_eidx(j,1),1)==p_vid(2)),  Map_vid_orient(j,1) = 1;  else    Map_vid_orient(j,1) = -1;  end
        if(Mesh.Edges(Map_vid_eidx(j,2),1)==p_vid(3)),  Map_vid_orient(j,2) = 1;  else    Map_vid_orient(j,2) = -1;  end
        if(Mesh.Edges(Map_vid_eidx(j,3),1)==p_vid(1)),  Map_vid_orient(j,3) = 1;  else    Map_vid_orient(j,3) = -1;  end
        
        % Compute element mapping
        p_bK= p_a1;
        p_BK = [p_a2-p_bK; ...
                    p_a3-p_bK];
        p_det_BK = abs(det(p_BK));
        p_inv_BK = inv(p_BK);
        p_TK = transpose(p_inv_BK); 
        
        point_hat = (point-p_bK)*p_inv_BK;
        Map_vid_values(j,:) = shap_W1F(point_hat);
        Map_vid_values(j,[1 2]) = Map_vid_values(j,[1 2])*p_TK*WM;
        Map_vid_values(j,[3 4]) = Map_vid_values(j,[3 4])*p_TK*WM;
        Map_vid_values(j,[5 6]) = Map_vid_values(j,[5 6])*p_TK*WM;
        
    end 
    
    % Determine the orientation
    if(Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else    p1 = -1;  end
    if(Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else    p3 = -1;  end
    
    % Add contributions to stiffness matrix
    % Aloc1 contribution shapefunction connected to first quadpoint
    Aloc1 = ...
        [(p1*N(1,[1,2])*reshape(Map_vid_values(1,:),2,3)).*Map_vid_orient(1,:); ...
         (p2*N(1,[3,4])*reshape(Map_vid_values(1,:),2,3)).*Map_vid_orient(1,:); ... 
         (p3*N(1,[5,6])*reshape(Map_vid_values(1,:),2,3)).*Map_vid_orient(1,:)];
    % Aloc2 contribution shapefunction connected to second quadpoint
    Aloc2= ...
        [(p1*N(2,[1,2])*reshape(Map_vid_values(2,:),2,3)).*Map_vid_orient(2,:); ...
         (p2*N(2,[3,4])*reshape(Map_vid_values(2,:),2,3)).*Map_vid_orient(2,:); ... 
         (p3*N(2,[5,6])*reshape(Map_vid_values(2,:),2,3)).*Map_vid_orient(2,:)];
    % Aloc3 contribution shapefunction connected to third quadpoint
    Aloc3 = ...
        [(p1*N(3,[1,2])*reshape(Map_vid_values(3,:),2,3)).*Map_vid_orient(3,:); ...
         (p2*N(3,[3,4])*reshape(Map_vid_values(3,:),2,3)).*Map_vid_orient(3,:); ... 
         (p3*N(3,[5,6])*reshape(Map_vid_values(3,:),2,3)).*Map_vid_orient(3,:)];    
    
    Aloc1 = Aloc1*det_BK/6;
    Aloc2 = Aloc2*det_BK/6;
    Aloc3 = Aloc3*det_BK/6;
    
    I(loc) = set_Rows(eidx,3);
    J(loc) = set_Cols(Map_vid_eidx(1,:),3);
    A(loc) = Aloc1(:);
    loc = loc+9;
    
    I(loc) = set_Rows(eidx,3);
    J(loc) = set_Cols(Map_vid_eidx(2,:),3);
    A(loc) = Aloc2(:);
    loc = loc+9;
    
    I(loc) = set_Rows(eidx,3);
    J(loc) = set_Cols(Map_vid_eidx(3,:),3);
    A(loc) = Aloc3(:);
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
