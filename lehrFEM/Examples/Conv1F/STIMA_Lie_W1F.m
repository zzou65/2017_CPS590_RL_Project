function Mloc = STIMA_Lie_W1F(Vertices,ElemInfo,V_HANDLE,QuadRule,varargin)
% STIMA_ContrRot computes element contribution of v x curl u term 
%
%   MLOC = STIMA_ContrRot(VERTICES ...) computes element contribution of v x curl u term matrix using 
%   Whitney 1-forms finite elements.
%
%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping
  Mloc=zeros(3,3);
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ]; % transpose of transformation matrix
  det_BK = abs(det(BK));     % twice the area of the triagle
  
  % Compute constant gradients of barycentric coordinate functions
  e1 = P3-P2;
  e2 = P1-P3;
  e3 = P2-P1;
  
  % Evaluate coefficient function at quadrature nodes
  Fval = -V_HANDLE([P1; P2; P3],ElemInfo,varargin{:});
  
   % Get barycentric coordinates of quadrature points 
  baryc= [1 0 0;  0 1 0; 0 0 1];
  
  % Compute constant gradients of barycentric coordinate functions
  g1 = [P2(2)-P3(2);P3(1)-P2(1)]/det_BK;
  g2 = [P3(2)-P1(2);P1(1)-P3(1)]/det_BK;
  g3 = [P1(2)-P2(2);P2(1)-P1(1)]/det_BK;
  
  % Evaluate basis functions at quadrature points
  % the rows of b(i) store the value of the the i-th
  % basis function at the quadrature points
  b1 = baryc(:,2)*g3'-baryc(:,3)*g2';
  b2 = baryc(:,3)*g1'-baryc(:,1)*g3';
  b3 = baryc(:,1)*g2'-baryc(:,2)*g1';
  
  r1=det( [e1; Fval(2,:) + Fval(3,:)]);
  r2=det( [e2; Fval(3,:) + Fval(1,:)]);
  r3=det( [e3; Fval(1,:) + Fval(2,:)]);
  if ( r1 >= 0 )
     Mloc(1,:) = r1/det_BK*[1 1 1]-...
         [Fval(3,:)*b1(3,:)'-Fval(2,:)*b1(2,:)', ...
          Fval(3,:)*b2(3,:)'-Fval(2,:)*b2(2,:)', ...
          Fval(3,:)*b3(3,:)'-Fval(2,:)*b3(2,:)'];
  end
  if ( r2 >= 0 )
      Mloc(2,:) = r2/det_BK*[1 1 1]-...
         [Fval(1,:)*b1(1,:)'-Fval(3,:)*b1(3,:)', ...
          Fval(1,:)*b2(1,:)'-Fval(3,:)*b2(3,:)', ...
          Fval(1,:)*b3(1,:)'-Fval(3,:)*b3(3,:)'];
  end
  if ( r3 >= 0 )
      Mloc(3,:) = r3/det_BK*[1 1 1]-...
         [Fval(2,:)*b1(2,:)'-Fval(1,:)*b1(1,:)', ...
          Fval(2,:)*b2(2,:)'-Fval(1,:)*b2(1,:)', ...
          Fval(2,:)*b3(2,:)'-Fval(1,:)*b3(1,:)'];
  end
  %Mloc=Mloc;
  return