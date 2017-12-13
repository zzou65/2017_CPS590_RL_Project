function Mloc = STIMA_ContrRot_UP(Vertices,ElemInfo,V_HANDLE,QuadRule,varargin)
% STIMA_ContrRot_UP computes element contribution of -v x curl u term in
% the sense of cochains, e.g.  mapping of edge DOFS to edge DOFS in
% evaluating  the intergralvalues \int_{edge_1} v x curl b_edge_j terms.
% Here we use the trace from the upwind direction.
%
%   MLOC = STIMA_ContrRot_UP(VERTICES ...) computes element contribution of v x curl u term matrix using 
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
  
  r1=det( [e1; Fval(2,:) + Fval(3,:)]);
  r2=det( [e2; Fval(3,:) + Fval(1,:)]);
  r3=det( [e3; Fval(1,:) + Fval(2,:)]);
  if ( r1 >= 0)
     Mloc(1,:) = r1/det_BK*[1 1 1];
  end
  if ( r2 >= 0)
      Mloc(2,:) = r2/det_BK*[1 1 1];
  end
  if ( r3 >= 0)
      Mloc(3,:) = r3/det_BK*[1 1 1];
  end
  %Mloc=Mloc;
  return