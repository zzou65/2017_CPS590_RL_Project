function Aloc = STIMA_LaplTensor_LFE(Vertices,ElemInfo,EHandle,QuadRule,varargin)
% STIMA_LaplTensor_LFE Element stiffness matrix for the Laplace with Tensor coefficiant equation.
%
%   ALOC = STIMA_LaplTensor(VERTICES,ELEMINFO,QUADRULE,FHANDLE) computes the 
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
%   EHANDLE is the function handle to the Tensor M_11 M_12 M_21 M_22
%   
%   ALOC = STIMA_LaplTensor_LFE(VERTICES,ELEMINFO,QUADRULE,FHANDLE,FPARAM) also
%   handles the additional variable length argument list FPARAM to the
%   function handle FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)[1 1 1 1];   
%   QuadRule = P7O6_rule();
%   Aloc = STIMA_LaplTensor([0 0; 1 0; 0 1],0,QuadRule,FHandle);
%
%   See also grad_shap_LFE.

%   Copyright 2005-2007 Patrick Meury & Kah Ling Sia & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPoints = size(QuadRule.w,1);
  
  % Preallocate memory
  
  Aloc = zeros(3,3);
  
  % Compute element mapping
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
  inv_BK = inv(BK);
  det_BK = abs(det(BK));
  
  TK = det_BK*transpose(inv_BK)*inv_BK;
  
  x = QuadRule.x*BK+ones(nPoints,1)*bK;
  
  % Compute element stiffness matrix
  
  FVal = EHandle(x,varargin{:});
  grad_N = grad_shap_LFE(QuadRule.x);
  
  grad_N_tr(:,1:2)=grad_N(:,1:2)*transpose(inv_BK);
  grad_N_tr(:,3:4)=grad_N(:,3:4)*transpose(inv_BK);
  grad_N_tr(:,5:6)=grad_N(:,5:6)*transpose(inv_BK);
  
  Aloc(1,1) = sum(QuadRule.w.*grad_N_tr(:,1).*FVal(:,1).*grad_N_tr(:,1))+...
              sum(QuadRule.w.*grad_N_tr(:,1).*FVal(:,2).*grad_N_tr(:,2))+...
              sum(QuadRule.w.*grad_N_tr(:,2).*FVal(:,3).*grad_N_tr(:,1))+...
              sum(QuadRule.w.*grad_N_tr(:,2).*FVal(:,4).*grad_N_tr(:,2));
  Aloc(1,2) = sum(QuadRule.w.*grad_N_tr(:,3).*FVal(:,1).*grad_N_tr(:,1))+...
              sum(QuadRule.w.*grad_N_tr(:,3).*FVal(:,2).*grad_N_tr(:,2))+...
              sum(QuadRule.w.*grad_N_tr(:,4).*FVal(:,3).*grad_N_tr(:,1))+...
              sum(QuadRule.w.*grad_N_tr(:,4).*FVal(:,4).*grad_N_tr(:,2));
  Aloc(1,3) = sum(QuadRule.w.*grad_N_tr(:,5).*FVal(:,1).*grad_N_tr(:,1))+...
              sum(QuadRule.w.*grad_N_tr(:,5).*FVal(:,2).*grad_N_tr(:,2))+...
              sum(QuadRule.w.*grad_N_tr(:,6).*FVal(:,3).*grad_N_tr(:,1))+...
              sum(QuadRule.w.*grad_N_tr(:,6).*FVal(:,4).*grad_N_tr(:,2));
  
  Aloc(2,1) = sum(QuadRule.w.*grad_N_tr(:,1).*FVal(:,1).*grad_N_tr(:,3))+...
              sum(QuadRule.w.*grad_N_tr(:,1).*FVal(:,2).*grad_N_tr(:,4))+...
              sum(QuadRule.w.*grad_N_tr(:,2).*FVal(:,3).*grad_N_tr(:,3))+...
              sum(QuadRule.w.*grad_N_tr(:,2).*FVal(:,4).*grad_N_tr(:,4));
  Aloc(2,2) = sum(QuadRule.w.*grad_N_tr(:,3).*FVal(:,1).*grad_N_tr(:,3))+...
              sum(QuadRule.w.*grad_N_tr(:,3).*FVal(:,2).*grad_N_tr(:,4))+...
              sum(QuadRule.w.*grad_N_tr(:,4).*FVal(:,3).*grad_N_tr(:,3))+...
              sum(QuadRule.w.*grad_N_tr(:,4).*FVal(:,4).*grad_N_tr(:,4));
  Aloc(2,3) = sum(QuadRule.w.*grad_N_tr(:,5).*FVal(:,1).*grad_N_tr(:,3))+...
              sum(QuadRule.w.*grad_N_tr(:,5).*FVal(:,2).*grad_N_tr(:,4))+...
              sum(QuadRule.w.*grad_N_tr(:,6).*FVal(:,3).*grad_N_tr(:,3))+...
              sum(QuadRule.w.*grad_N_tr(:,6).*FVal(:,4).*grad_N_tr(:,4));
  
  Aloc(3,1) = sum(QuadRule.w.*grad_N_tr(:,1).*FVal(:,1).*grad_N_tr(:,5))+...
              sum(QuadRule.w.*grad_N_tr(:,1).*FVal(:,2).*grad_N_tr(:,6))+...
              sum(QuadRule.w.*grad_N_tr(:,2).*FVal(:,3).*grad_N_tr(:,5))+...
              sum(QuadRule.w.*grad_N_tr(:,2).*FVal(:,4).*grad_N_tr(:,6));
  Aloc(3,2) = sum(QuadRule.w.*grad_N_tr(:,3).*FVal(:,1).*grad_N_tr(:,5))+...
              sum(QuadRule.w.*grad_N_tr(:,3).*FVal(:,2).*grad_N_tr(:,6))+...
              sum(QuadRule.w.*grad_N_tr(:,4).*FVal(:,3).*grad_N_tr(:,5))+...
              sum(QuadRule.w.*grad_N_tr(:,4).*FVal(:,4).*grad_N_tr(:,6));
  Aloc(3,3) = sum(QuadRule.w.*grad_N_tr(:,5).*FVal(:,1).*grad_N_tr(:,5))+...
              sum(QuadRule.w.*grad_N_tr(:,5).*FVal(:,2).*grad_N_tr(:,6))+...
              sum(QuadRule.w.*grad_N_tr(:,6).*FVal(:,3).*grad_N_tr(:,5))+...
              sum(QuadRule.w.*grad_N_tr(:,6).*FVal(:,4).*grad_N_tr(:,6));
          
 Aloc=det_BK*Aloc;
  
return