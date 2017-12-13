function varargout= assemMat_SemiLagQuad_W1F_strang12nd(Mesh,velocity,Jac,h,varargin)
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

options = odeset('MaxStep',h);
nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

% Preallocate memory
nElements = size(Mesh.Elements,1);

% Preallocate memory

I = zeros(3*36*nElements,1);
J = zeros(3*36*nElements,1);
A = zeros(3*36*nElements,1);

% Check for element flags
if (isfield(Mesh,'ElemFlag')), flags = Mesh.ElemFlag;
else flags = zeros(nElements,1); end

% Assemble element contributions
loc = 1:36;
for i = 1:nElements
   
    Aloc1=zeros(6,6);
    Aloc2=zeros(6,6);
    Aloc3=zeros(6,6);

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

    % loop over all three quadraturepoints
    Map_vid_eidx = zeros(3,3);
    Map_vid_values = zeros(3,12);
    Map_vid_orient = ones(3,6);
    
    % quadraturepoints and weights
    x =1/6*[ (a1+a2+4*a3); (4*a1+a2+a3); (a1+4*a2+a3)];
    w =1/3*[1 1 1];
    
    N = shap_W1F2nd([1/6 4/6; 1/6 1/6; 4/6 1/6]);
    N(:,[1 2]) = N(:,[1 2])*TK;
    N(:,[3 4]) = N(:,[3 4])*TK;
    N(:,[5 6]) = N(:,[5 6])*TK;
    N(:,[7 8]) = N(:,[7 8])*TK;
    N(:,[9 10]) = N(:,[9 10])*TK;
    N(:,[11 12]) = N(:,[11 12])*TK;
        
    for j=1:3
        % traced quadrature point and Jacobian
        %rhs = @(t,x)[velocity([x(1) x(2)])'; ...
        %    sum(Jac([x(1) x(2)]).*[x(3) 0 x(4) 0],2); ...
        %    sum(Jac([x(1) x(2)]).*[0 x(3) 0 x(4)],2); ...
        %    sum(Jac([x(1) x(2)]).*[x(5) 0 x(6) 0],2); ...
        %    sum(Jac([x(1) x(2)]).*[0 x(5) 0 x(6)],2)];
        %y0 = [x(j,:)'; 1 ; 0; 0; 1];
        %[t,y] = ode45(rhs,[0,h],y0,options); 
        %y=y0'+h*rhs(0,y0)'; 
        %dir=y(end,[1 2])-x(j,:);
        dir = h*velocity(x(j,:);
        traced_qpoint=trace_ipoint(i,x(j,:),dir,Mesh);
        
        WM=[1 0 0 1] + h*Jac(x(j,:)); 

        point=traced_qpoint([1,2]);
        p_Element=traced_qpoint(3);
        
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
        Map_vid_values(j,:) = shap_W1F2nd(point_hat);
        Map_vid_values(j,[1 2]) = Map_vid_values(j,[1 2])*p_TK*WM;
        Map_vid_values(j,[3 4]) = Map_vid_values(j,[3 4])*p_TK*WM;
        Map_vid_values(j,[5 6]) = Map_vid_values(j,[5 6])*p_TK*WM;
        Map_vid_values(j,[7 8]) = Map_vid_values(j,[7  8])*p_TK*WM;
        Map_vid_values(j,[9 10]) = Map_vid_values(j,[9 10])*p_TK*WM;
        Map_vid_values(j,[11 12]) = Map_vid_values(j,[11 12])*p_TK*WM;
    end 
    % Determine the orientation

    if(Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else    p1 = -1;  end
    if(Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else    p3 = -1;  end
    
    % Add contributions to stiffness matrix
    % Aloc1 contribution shapefunction connected to first quadpoint
    Aloc1 = ...
        [(p1*N(1,[1,2])*reshape(Map_vid_values(1,:),2,6)).*Map_vid_orient(1,:); ...
         (p2*N(1,[3,4])*reshape(Map_vid_values(1,:),2,6)).*Map_vid_orient(1,:); ... 
         (p3*N(1,[5,6])*reshape(Map_vid_values(1,:),2,6)).*Map_vid_orient(1,:); ... 
         (N(1,[7,8])*reshape(Map_vid_values(1,:),2,6)).*Map_vid_orient(1,:); ... 
         (N(1,[9,10])*reshape(Map_vid_values(1,:),2,6)).*Map_vid_orient(1,:); ... 
         (N(1,[11,12])*reshape(Map_vid_values(1,:),2,6)).*Map_vid_orient(1,:)
         ];
    % Aloc2 contribution shapefunction connected to second quadpoint
    Aloc2= ...
        [(p1*N(2,[1,2])*reshape(Map_vid_values(2,:),2,6)).*Map_vid_orient(2,:); ...
         (p2*N(2,[3,4])*reshape(Map_vid_values(2,:),2,6)).*Map_vid_orient(2,:); ... 
         (p3*N(2,[5,6])*reshape(Map_vid_values(2,:),2,6)).*Map_vid_orient(2,:); ... 
         (N(2,[7,8])*reshape(Map_vid_values(2,:),2,6)).*Map_vid_orient(2,:); ... 
         (N(2,[9,10])*reshape(Map_vid_values(2,:),2,6)).*Map_vid_orient(2,:); ... 
         (N(2,[11,12])*reshape(Map_vid_values(2,:),2,6)).*Map_vid_orient(2,:)];
    % Aloc3 contribution shapefunction connected to third quadpoint
    Aloc3 = ...
        [(p1*N(3,[1,2])*reshape(Map_vid_values(3,:),2,6)).*Map_vid_orient(3,:); ...
         (p2*N(3,[3,4])*reshape(Map_vid_values(3,:),2,6)).*Map_vid_orient(3,:); ... 
         (p3*N(3,[5,6])*reshape(Map_vid_values(3,:),2,6)).*Map_vid_orient(3,:); ... 
         (N(3,[7,8])*reshape(Map_vid_values(3,:),2,6)).*Map_vid_orient(3,:); ... 
         (N(3,[9,10])*reshape(Map_vid_values(3,:),2,6)).*Map_vid_orient(3,:); ... 
         (N(3,[11,12])*reshape(Map_vid_values(3,:),2,6)).*Map_vid_orient(3,:)];    
    
    Aloc1 = w(1)*Aloc1*det_BK/2;
    Aloc2 = w(2)*Aloc2*det_BK/2;
    Aloc3 = w(3)*Aloc3*det_BK/2;
    
    I(loc) = set_Rows([eidx nEdges+eidx],6);
    J(loc) = set_Cols([Map_vid_eidx(1,:) nEdges+Map_vid_eidx(1,:)],6);
    A(loc) = Aloc1(:);
    loc = loc+36;
    
    I(loc) = set_Rows([eidx nEdges+eidx],6);
    J(loc) = set_Cols([Map_vid_eidx(2,:) nEdges+Map_vid_eidx(2,:)],6);
    A(loc) = Aloc2(:);
    loc = loc+36;
    
    I(loc) = set_Rows([eidx nEdges+eidx],6);
    J(loc) = set_Cols([Map_vid_eidx(3,:) nEdges+Map_vid_eidx(3,:)],6);
    A(loc) = Aloc3(:);
    loc = loc+36;
    
end

% Assign output arguments

if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I,J,A,2*nEdges,2*nEdges);
end

function f=rhs(x,t)
    f = velocity(x);
