function varargout = assemMat_ContrRot_UPQuad(Mesh,V_HANDLE,varargin)
% consistent scheme for -v x curl u based on upwind quadrature using P3O3() 

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  [EdgeContr C]=assemData_ContrRot_UpQuad(Mesh, V_HANDLE);

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  I = zeros(3*9*nElements,1);
  J = zeros(3*9*nElements,1);
  A = zeros(3*9*nElements,1);
  
  QuadRule=P3O3();
  % Assemble element contributions
  
  loc = 1:9;
  for i = 1:nElements
    
    % Extract vertices of current element
    
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
    
    % Extract global edge numbers
    
    eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];
    
    % Compute element contributions

    P1 = Vertices(1,:);
    P2 = Vertices(2,:);
    P3 = Vertices(3,:);

    BK = [ P2 - P1 ; P3 - P1 ];   % transpose of transformation matrix
    det_BK = abs(det(BK));        % twice the area of the triagle

    % Compute constant gradients of barycentric coordinate functions
    g1 = [P2(2)-P3(2);P3(1)-P2(1)]/det_BK;
    g2 = [P3(2)-P1(2);P1(1)-P3(1)]/det_BK;
    g3 = [P1(2)-P2(2);P2(1)-P1(1)]/det_BK;

    % Get barycentric coordinates of quadrature points
    nPoints = size(QuadRule.w,1);
    baryc= [1-sum(QuadRule.x,2),QuadRule.x];

    % Quadrature points in actual element
    % stored as rows of a matrix
    x = QuadRule.x*BK + ones(nPoints,1)*P1;

    % Evaluate coefficient function at quadrature nodes
    Fval = V_HANDLE(x);

    %Rotation
    Fval =[Fval(:,2) -Fval(:,1)];

    % Evaluate basis functions at quadrature points
    % the rows of b(i) store the value of the the i-th
    % basis function at the quadrature points
    b1 = baryc(:,2)*g3'-baryc(:,3)*g2';
    b2 = baryc(:,3)*g1'-baryc(:,1)*g3';
    b3 = baryc(:,1)*g2'-baryc(:,2)*g1';

    % Determine the orientation
    
    if(Mesh.Edges(eidx(1),1)==vidx(2)),  p1 = 1;  else    p1 = -1;  end
    if(Mesh.Edges(eidx(2),1)==vidx(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vidx(1)),  p3 = 1;  else    p3 = -1;  end
        
    % Compute local mass matrix
    weights = det_BK/6;
    % first quad point
    Aloc(1,1) = [p1*weights.*sum(Fval(2,:).*b1(2,:),2).*EdgeContr(eidx(1),1)];
    Aloc(1,2) = [p1*weights.*sum(Fval(2,:).*b1(2,:),2).*EdgeContr(eidx(1),2)];
    Aloc(1,3) = [p1*weights.*sum(Fval(2,:).*b1(2,:),2).*EdgeContr(eidx(1),3)];

    Aloc(2,1) = [p2*weights.*sum(Fval(2,:).*b2(2,:),2).*EdgeContr(eidx(1),1)];
    Aloc(2,2) = [p2*weights.*sum(Fval(2,:).*b2(2,:),2).*EdgeContr(eidx(1),2)];
    Aloc(2,3) = [p2*weights.*sum(Fval(2,:).*b2(2,:),2).*EdgeContr(eidx(1),3)];

    Aloc(3,1) = [p3*weights.*sum(Fval(2,:).*b3(2,:),2).*EdgeContr(eidx(1),1)];
    Aloc(3,2) = [p3*weights.*sum(Fval(2,:).*b3(2,:),2).*EdgeContr(eidx(1),2)];
    Aloc(3,3) = [p3*weights.*sum(Fval(2,:).*b3(2,:),2).*EdgeContr(eidx(1),3)];
    
    % Add contributions to stiffness matrix
    if C(eidx(1),1)~=0
        I(loc) = set_Rows(eidx,3);
        J(loc) = set_Cols(C(eidx(1),:),3);
        A(loc) = -Aloc(:);
        loc = loc+9;
    end
    
     % second quad point
    Aloc(1,1) = [p1*weights.*sum(Fval(1,:).*b1(1,:),2).*EdgeContr(eidx(2),1)];
    Aloc(1,2) = [p1*weights.*sum(Fval(1,:).*b1(1,:),2).*EdgeContr(eidx(2),2)];
    Aloc(1,3) = [p1*weights.*sum(Fval(1,:).*b1(1,:),2).*EdgeContr(eidx(2),3)];

    Aloc(2,1) = [p2*weights.*sum(Fval(1,:).*b2(1,:),2).*EdgeContr(eidx(2),1)];
    Aloc(2,2) = [p2*weights.*sum(Fval(1,:).*b2(1,:),2).*EdgeContr(eidx(2),2)];
    Aloc(2,3) = [p2*weights.*sum(Fval(1,:).*b2(1,:),2).*EdgeContr(eidx(2),3)];

    Aloc(3,1) = [p3*weights.*sum(Fval(1,:).*b3(1,:),2).*EdgeContr(eidx(2),1)];
    Aloc(3,2) = [p3*weights.*sum(Fval(1,:).*b3(1,:),2).*EdgeContr(eidx(2),2)];
    Aloc(3,3) = [p3*weights.*sum(Fval(1,:).*b3(1,:),2).*EdgeContr(eidx(2),3)];
    
    % Add contributions to stiffness matrix
    if C(eidx(2),1)~=0
        I(loc) = set_Rows(eidx,3);
        J(loc) = set_Cols(C(eidx(2),:),3);
        A(loc) = -Aloc(:);
        loc = loc+9;
    end
    
   % third quad point
    Aloc(1,1) = [p1*weights.*sum(Fval(3,:).*b1(3,:),2).*EdgeContr(eidx(3),1)];
    Aloc(1,2) = [p1*weights.*sum(Fval(3,:).*b1(3,:),2).*EdgeContr(eidx(3),2)];
    Aloc(1,3) = [p1*weights.*sum(Fval(3,:).*b1(3,:),2).*EdgeContr(eidx(3),3)];

    Aloc(2,1) = [p2*weights.*sum(Fval(3,:).*b2(3,:),2).*EdgeContr(eidx(3),1)];
    Aloc(2,2) = [p2*weights.*sum(Fval(3,:).*b2(3,:),2).*EdgeContr(eidx(3),2)];
    Aloc(2,3) = [p2*weights.*sum(Fval(3,:).*b2(3,:),2).*EdgeContr(eidx(3),3)];

    Aloc(3,1) = [p3*weights.*sum(Fval(3,:).*b3(3,:),2).*EdgeContr(eidx(3),1)];
    Aloc(3,2) = [p3*weights.*sum(Fval(3,:).*b3(3,:),2).*EdgeContr(eidx(3),2)];
    Aloc(3,3) = [p3*weights.*sum(Fval(3,:).*b3(3,:),2).*EdgeContr(eidx(3),3)];
    
    % Add contributions to stiffness matrix
    if C(eidx(3),1)~=0
        I(loc) = set_Rows(eidx,3);
        J(loc) = set_Cols(C(eidx(3),:),3);
        A(loc) = -Aloc(:);
        loc = loc+9;
    end
    
  end
  
  % Assign output arguments
  I=I(1:loc-9);
  J=J(1:loc-9);
  A=A(1:loc-9);
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A);      
  end
  
return