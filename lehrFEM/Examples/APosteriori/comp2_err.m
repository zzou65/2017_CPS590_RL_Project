function err = comp2_err(Mesh,u,QuadRule,FHandle,varargin)
% computes the integral of u-FHandle; based on L1Err_LFE but without the
% absolute value inside the integral; used for the setting present in
% comp_2.m

%   Copyright 2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  
  % Precompute shape functions
  
  N = shap_LFE(QuadRule.x);
    
  % Compute discretization error

  err = 0;
  for i = 1:nElements
       
    % Extract vertex numbers
    
    vidx = Mesh.Elements(i,:);
      
    % Compute element mapping  
        
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));

    % Transform quadrature points
      
    x = QuadRule.x*BK+ones(nPts,1)*bK;
      
    % Evaluate solutions
      
    u_EX = FHandle(x,varargin{:});
    u_FE = u(vidx(1))*N(:,1) + u(vidx(2))*N(:,2) + u(vidx(3))*N(:,3);
      
    % Compute error on current element
    
    err = err+sum(QuadRule.w.*(u_EX-u_FE))*det_BK;
      
  end
  
return