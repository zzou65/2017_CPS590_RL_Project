function Aloc = STIMA_Curl_W1F2nd(Vertices,ElemInfo,MU_HANDLE,QuadRule,varargin)
% STIMA_CURL_W1F element stiffness matrix for curl*curl-operator in 2D
% in the case of Galerkin discretization by means of edge elements
%
%   ALOC = STIMA_CURL_W1F(VERTICES,ELEMINFO,MU_HANDLE,QUADRULE) computes the
%   curl*\mu*curl element stiffness matrix using Whitney 1-forms finite elements.
%   The function \mu can be passed through the MU_HANDLE argument
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%   
%   ElemInfo (not used)
%
%   MU_HANDLE handle to a functions expecting a matrix whose rows
%   represent position arguments. Return value must be a vector
%   (variable arguments will be passed to this function)
%   
%  QuadRule is a quadrature rule on the reference element
%
%   Example:
%
%   Aloc = STIMA_Curl_W1F(Vertices,ElemInfo,MU_HANDLE,QuadRule);

%   Copyright 2005-2006Patrick Meury & Mengyu Wang & Ralf Hiptmair
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constant
  
  nPoints = size(QuadRule.w,1);
  
  % Compute element mapping

  bK = Vertices(1,:); % row vector !
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK]; % Transpose of trafo matrix !
  det_BK = abs(det(BK)); % twice the area of the triangle
  
  % Quadrature points in actual element
  % stored as rows of a matrix
  x = QuadRule.x*BK + ones(nPoints,1)*bK;
  
  % Compute function value
  
  Fval = MU_HANDLE(x,ElemInfo,varargin{:});
  
  % Compute local curl-curl-matrix
  % Use that the curl of an edge element function is constant
  % and equals 1/area of triangle
  
  Aloc = 4/det_BK*sum(QuadRule.w.*Fval)*[ones(3,3) zeros(3,3); zeros(3,3) zeros(3,3)] ;
return
