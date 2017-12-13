function Mloc =STIMA_ContrGrad_LFE(Vertices, dummy ,VHandle)
% MASS_LFE Element mass matrix.
%
%   Mloc = STIMA_ContrGrad(VERTICES,dummy, VHandle) computes the element mass matrix using linear
%   Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   Example:
%
%   Mloc = STIMA_ContrGrad_LFE(Vertices,0, Vhandle);

%   Copyright 2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Compute element mapping

  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
  det_BK = abs(det(BK));

% Compute local mass matrix
  
  v1=VHandle(Vertices(1,:));
  v2=VHandle(Vertices(2,:));
  v3=VHandle(Vertices(3,:));
  
% Compute gradient of barycentric coordinate functions
  
  grad_Lambda1=1/det_BK*[Vertices(2,2)-Vertices(3,2) -Vertices(2,1)+Vertices(3,1)]; 
  grad_Lambda2=1/det_BK*[Vertices(3,2)-Vertices(1,2) -Vertices(3,1)+Vertices(1,1)];
  grad_Lambda3=1/det_BK*[Vertices(1,2)-Vertices(2,2) -Vertices(1,1)+Vertices(2,1)];
  
% Compute components of the local element mass matrix using linear
% Lagrangian finite elements
  
  Mloc(1,1)=det_BK/6*grad_Lambda1*v1';
  Mloc(1,2)=det_BK/6*grad_Lambda1*v2';
  Mloc(1,3)=det_BK/6*grad_Lambda1*v3';
  
  Mloc(2,1)=det_BK/6*grad_Lambda2*v1';
  Mloc(2,2)=det_BK/6*grad_Lambda2*v2';
  Mloc(2,3)=det_BK/6*grad_Lambda2*v3';
  
  Mloc(3,1)=det_BK/6*grad_Lambda3*v1';
  Mloc(3,2)=det_BK/6*grad_Lambda3*v2';
  Mloc(3,3)=det_BK/6*grad_Lambda3*v3';
  
 Mloc=Mloc';
   
return
