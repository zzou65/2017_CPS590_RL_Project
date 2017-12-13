function varargout = assemMat_Contr_1FVertUP(Mesh,V_HANDLE,varargin)
%  assemMatContr_1FVerUP assembles Contraction based on upwind quadrature
%  at vertices

%   Copyright 2005-2009 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  D = assemDataVert_Contr1f(Mesh, V_HANDLE);

  % Initialize constants
  nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);
  % Preallocate memory
  
  I = zeros(6*nElements,1);
  J = zeros(6*nElements,1);
  A = zeros(6*nElements,1);
  
  QuadRule=P3O2();
  % Assemble element contributions
  
  loc = 0;
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

    % Quadrature points in actual element
    % stored as rows of a matrix
    x = QuadRule.x*BK + ones(3,1)*P1;

    % Evaluate coefficient function at quadrature nodes
    Fval = V_HANDLE(x);

    %Rotation
    Fval =[Fval(:,2) -Fval(:,1)];

    % Determine the orientation
    
    if(Mesh.Edges(eidx(1),1)==vidx(2)),  p1 = 1;  else    p1 = -1;  end
    if(Mesh.Edges(eidx(2),1)==vidx(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vidx(1)),  p3 = 1;  else    p3 = -1;  end
        
    % Compute local mass matrix
     
    if D(vidx(1),1)~=0
        %second edge
        A(loc+2) = -p2*1/4*D(vidx(1),2)*Fval(1,:)*(P1-P3)';
        I(loc+2) = eidx(2);
        J(loc+2) = D(vidx(1),1);
        
        % third edge
        A(loc+1) = -p3*1/4*D(vidx(1),2)*Fval(1,:)*(P2-P1)';
        I(loc+1) = eidx(3);
        J(loc+1) = D(vidx(1),1);
        
        loc=loc+2;
    end
    if D(vidx(2),1)~=0
        % first edge
        A(loc+1) = -p1*1/4*D(vidx(2),2)*Fval(2,:)*(P3-P2)';
        I(loc+1) = eidx(1);
        J(loc+1) = D(vidx(2),1);

        % third edge
        A(loc+2) =- p3*1/4*D(vidx(2),2)*Fval(2,:)*(P2-P1)';
        I(loc+2) = eidx(3);
        J(loc+2) = D(vidx(2),1);
        
        loc=loc+2;
    end
    
    if D(vidx(3),1)~=0
        % first edge 
        A(loc+2) = -p1*1/4*D(vidx(3),2)*Fval(3,:)*(P3-P2)';
        I(loc+2) = eidx(1);
        J(loc+2) = D(vidx(3),1);

        % second edge
        A(loc+1) = -p2*1/4*D(vidx(3),2)*Fval(3,:)*(P1-P3)';
        I(loc+1) = eidx(2);
        J(loc+1) = D(vidx(3),1);

        loc=loc+2;
    end
  end
  I=I(1:loc-2);
  J=J(1:loc-2);
  A=A(1:loc-2);
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A,nEdges,nElements);      
  end
  
return