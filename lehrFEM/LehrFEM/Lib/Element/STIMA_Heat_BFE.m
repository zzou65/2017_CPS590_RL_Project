function Aloc = STIMA_Heat_BFE(Vertices,ElemInfo,QuadRule,FHandle,varargin)
% STIMA_HEAT_BFE Element stiffness matrix for the Laplacian.
%
%   ALOC = STIMA_HEAT_BFE(VERTICES,ELEMINFO,QUADRULE,FHANDLE) computes the 
%   element stiffness matrix for the data given by function handle FHANDLE. 
%
%   VERTICES is a 4-by-2 matrix specifying the vertices of the current
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
%   ALOC = STIMA_HEAT_BFE(VERTICES,ELEMINFO,QUADRULE,FHANDLE,FPARAM) also
%   handles the additional variable length argument list FPARAM to the
%   function handle FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)1;   
%   QuadRule = TProd(gauleg(0,1,NGAUSS));
%   Aloc = STIMA_Heat_BFE([0 0; 1 0; 1 1; 0 1],0,QuadRule,FHandle);
%
%   See also shap_BFE, grad_shap_BFE.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nGauss = size(QuadRule.w,1);
  
  % Preallocate memory
  
  Aloc = zeros(4,4);
  DPhi_K = zeros(2,2);
  
  % Compute values of shape functions
  
  shap_N = shap_BFE(QuadRule.x);
  grad_N = grad_shap_BFE(QuadRule.x);
  
  % Extract vertices
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  P4 = Vertices(4,:);
  
  % Compute function values
  
  x = shap_N(:,1)*P1+shap_N(:,2)*P2+shap_N(:,3)*P3+shap_N(:,4)*P4;
  FVal = FHandle(x,ElemInfo,varargin{:});
      
  for i = 1:nGauss

    % Compute element map  
    
    DPhi_K(1,:) = P1(1)*grad_N(i,1:2)+P2(1)*grad_N(i,3:4) + ...
                  P3(1)*grad_N(i,5:6)+P4(1)*grad_N(i,7:8);
    DPhi_K(2,:) = P1(2)*grad_N(i,1:2)+P2(2)*grad_N(i,3:4) + ...
                  P3(2)*grad_N(i,5:6)+P4(2)*grad_N(i,7:8);
    inv_DPhi_K = inv(DPhi_K);
    TK = inv_DPhi_K*transpose(inv_DPhi_K)*abs(det(DPhi_K));
    
    % Compute entries of stiffness matrix
    
    Aloc(1,1) = Aloc(1,1) + FVal(i)*QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,1:2));
    Aloc(1,2) = Aloc(1,2) + FVal(i)*QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,3:4));
    Aloc(1,3) = Aloc(1,3) + FVal(i)*QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,5:6));
    Aloc(1,4) = Aloc(1,4) + FVal(i)*QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,7:8));
    
    Aloc(2,2) = Aloc(2,2) + FVal(i)*QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,3:4));
    Aloc(2,3) = Aloc(2,3) + FVal(i)*QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,5:6));
    Aloc(2,4) = Aloc(2,4) + FVal(i)*QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,7:8));
    
    Aloc(3,3) = Aloc(3,3) + FVal(i)*QuadRule.w(i)*grad_N(i,5:6)*TK*transpose(grad_N(i,5:6));
    Aloc(3,4) = Aloc(3,4) + FVal(i)*QuadRule.w(i)*grad_N(i,5:6)*TK*transpose(grad_N(i,7:8));
    
    Aloc(4,4) = Aloc(4,4) + FVal(i)*QuadRule.w(i)*grad_N(i,7:8)*TK*transpose(grad_N(i,7:8));
    
  end
  
  % Fill in lower triangular part
  
  tri = triu(Aloc);
  Aloc = tri+tril(tri',-1);
  
return