function Aloc = STIMA_WReg_W1Fa(Vertices,ElemInfo,varargin)
% STIMA_WREG_W1F Element stiffness matrix for the W1F finite element.
%
%   ALOC = STIMA_WREG_W1F(VERTICES,ELEMINFO) computes the element stiffness 
%   matrix for the data given by function handle FHANDLE. 
%
%   VERTICES is a 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   Example:
%
%   Aloc = STIMA_WReg_W1F([0 0; 1 0; 0 1],0);
%
%   See also grad_shap_LFE.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Aloc = zeros(3,3);
  
  % Compute element mapping
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  bK = P1;
  BK = [P2-bK;P3-bK];
  inv_BK = inv(BK);
  det_BK = abs(det(BK));
  TK = transpose(inv_BK);
  
  phi_1 = 1/6*[2*P1(2)-P2(2)-P3(2) -2*P1(1)+P2(1)+P3(1)];
  phi_2 = 1/6*[2*P2(2)-P3(2)-P1(2) -2*P2(1)+P3(1)+P1(1)];
  phi_3 = 1/6*[2*P3(2)-P1(2)-P2(2) -2*P3(1)+P1(1)+P2(1)];
  
  % Compute element stiffness matrix

  Aloc(1,1) = dot([-1 -1]*TK,phi_1);
  Aloc(2,1) = dot([-1 -1]*TK,phi_2);
  Aloc(3,1) = dot([-1 -1]*TK,phi_3);
  Aloc(1,2) = dot([1 0]*TK,phi_1);
  Aloc(2,2) = dot([1 0]*TK,phi_2);
  Aloc(3,2) = dot([1 0]*TK,phi_3);
  Aloc(1,3) = dot([0,1]*TK,phi_1);
  Aloc(2,3) = dot([0,1]*TK,phi_2);
  Aloc(3,3) = dot([0 1]*TK,phi_3);
  
return