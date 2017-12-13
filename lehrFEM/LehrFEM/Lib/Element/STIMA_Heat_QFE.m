function Aloc = STIMA_Heat_QFE(Vertices,ElemInfo,QuadRule,FHandle,varargin)
% STIMA_HEAT_QFE Element stiffness matrix for the stationary heat equation.
%
%   ALOC = STIMA_HEAT_QFE(VERTICES,ELEMINFO,QUADRULE,FHANDLE) computes the 
%   element stiffness matrix for the data given by function handle FHANDLE. 
%
%   VERTICES is a 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%    
%   FHANDLE is the function handle to the heat conductivity of the
%   stationary heat conduction equation.
%   
%   ALOC = STIMA_HEAT_QFE(VERTICES,ELEMINFO,QUADRULE,FHANDLE,FPARAM) also
%   handles the additional variable length argument list FPARAM to the
%   function handle FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)1;   
%   QuadRule = P7O6_rule();
%   Aloc = STIMA_Heat_QFE([0 0; 1 0; 0 1],0,QuadRule,FHandle);
%
%   See also GRAD_SHAP_QFE.m

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  
  % Initialize constants
  
  nPoints = size(QuadRule.w,1);
  
  % Preallocate memory
  
  Aloc = zeros(6,6);
  
  % Compute element mapping
  
  BK = [Vertices(2,:)-Vertices(1,:); Vertices(3,:)-Vertices(1,:)];
  bK = Vertices(1,:);
  inv_BK = inv(BK);
  det_BK = abs(det(BK));
  
  TK = det_BK*transpose(inv_BK)*inv_BK;
  
  x = QuadRule.x*BK+ones(nPoints,1)*bK;
  
  % Compute element stiffness matrix
  
  FVal = FHandle(x,ElemInfo,varargin{:});
  grad_N = grad_shap_QFE(QuadRule.x);
  
  Aloc(1,1) = sum(QuadRule.w.*FVal.*sum((grad_N(:,1:2)).*(grad_N(:,1:2)*TK),2));
  Aloc(1,2) = sum(QuadRule.w.*FVal.*sum((grad_N(:,1:2)).*(grad_N(:,3:4)*TK),2));
  Aloc(1,3) = sum(QuadRule.w.*FVal.*sum((grad_N(:,1:2)).*(grad_N(:,5:6)*TK),2));
  Aloc(1,4) = sum(QuadRule.w.*FVal.*sum((grad_N(:,1:2)).*(grad_N(:,7:8)*TK),2));
  Aloc(1,5) = sum(QuadRule.w.*FVal.*sum((grad_N(:,1:2)).*(grad_N(:,9:10)*TK),2));
  Aloc(1,6) = sum(QuadRule.w.*FVal.*sum((grad_N(:,1:2)).*(grad_N(:,11:12)*TK),2));
  
  Aloc(2,2) = sum(QuadRule.w.*FVal.*sum((grad_N(:,3:4)).*(grad_N(:,3:4)*TK),2));
  Aloc(2,3) = sum(QuadRule.w.*FVal.*sum((grad_N(:,3:4)).*(grad_N(:,5:6)*TK),2));
  Aloc(2,4) = sum(QuadRule.w.*FVal.*sum((grad_N(:,3:4)).*(grad_N(:,7:8)*TK),2));
  Aloc(2,5) = sum(QuadRule.w.*FVal.*sum((grad_N(:,3:4)).*(grad_N(:,9:10)*TK),2));
  Aloc(2,6) = sum(QuadRule.w.*FVal.*sum((grad_N(:,3:4)).*(grad_N(:,11:12)*TK),2));
  
  Aloc(3,3) = sum(QuadRule.w.*FVal.*sum((grad_N(:,5:6)).*(grad_N(:,5:6)*TK),2));
  Aloc(3,4) = sum(QuadRule.w.*FVal.*sum((grad_N(:,5:6)).*(grad_N(:,7:8)*TK),2));
  Aloc(3,5) = sum(QuadRule.w.*FVal.*sum((grad_N(:,5:6)).*(grad_N(:,9:10)*TK),2));
  Aloc(3,6) = sum(QuadRule.w.*FVal.*sum((grad_N(:,5:6)).*(grad_N(:,11:12)*TK),2));
  
  Aloc(4,4) = sum(QuadRule.w.*FVal.*sum((grad_N(:,7:8)).*(grad_N(:,7:8)*TK),2));
  Aloc(4,5) = sum(QuadRule.w.*FVal.*sum((grad_N(:,7:8)).*(grad_N(:,9:10)*TK),2));
  Aloc(4,6) = sum(QuadRule.w.*FVal.*sum((grad_N(:,7:8)).*(grad_N(:,11:12)*TK),2));
  
  Aloc(5,5) = sum(QuadRule.w.*FVal.*sum((grad_N(:,9:10)).*(grad_N(:,9:10)*TK),2));
  Aloc(5,6) = sum(QuadRule.w.*FVal.*sum((grad_N(:,9:10)).*(grad_N(:,11:12)*TK),2));
  
  Aloc(6,6) = sum(QuadRule.w.*FVal.*sum((grad_N(:,11:12)).*(grad_N(:,11:12)*TK),2));
  
     
  % Fill in lower triangular part
  
  tri = triu(Aloc);
  Aloc = tri+tril(tri',-1);
  
return