function u = mySlopeLim_lr(Coordinates,v,p,Inn_shap, shap, QuadRule,m,varargin)
% less restictive slope limiter according to Cockburn & Shu, and a
% modification, that remains polynomial degree after restriction
%
%   Coordinates: n+1 Meshpoints
%             u:
%             p: polynomial degree in each of the n cells
%      Inn_shap: basisfunctions evaluated at cell boundary
%    flowhandle: flux, numericalflux 
%      Quadrule: 1D quadrature rule 
%             m: parameter for modified minmod
%   Example:
%   
%    
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
 % Initialize constants
  
  
  if (min(p)==0);
      u=v;
  else
  % Initialize constants
  NCoordinates=size(Coordinates,2);
  NElements=NCoordinates-1;
  Dofs=sum(p+1);
  lshap=Inn_shap(1,:);
  rshap=Inn_shap(2,:);
  % Alocate memory
  u=zeros(Dofs,1);
  
  % left cell average
  vm_l=v(sum(p(1:NElements-1)+1)+1);
  
  % middle cell average
  vm_m=v(1);
  
  offset=0;                    
  % loop over cells
  for i=1:(NElements-1)
     
      % Cellsize and dofs of middle cell
      dofs=p(i)+1;
      h=Coordinates(i+1)-Coordinates(i); 
      
      % right cell average
      vm_r=v(offset+dofs+1);
      
      % limiter calculation
      
      % left values, exact 
      v_pos=lshap(1:dofs)*v(offset+(1:dofs));
      
      % right values exact 
      v_neg=rshap(1:dofs)*v(offset+(1:dofs));
      
      u(offset+(1:dofs))=v(offset+(1:dofs));
      u(offset+2)=minmod(v(offset+2),...
             (vm_m-vm_l)-(v_pos-vm_m-v(offset+2)),...
             (vm_r-vm_m)-(v_neg-vm_m-v(offset+2)),m*h^2);
    
      offset=offset+dofs;
      vm_l=vm_m;
      vm_m=vm_r;
      
  end
  dofs=p(NElements)+1;
  h=Coordinates(NCoordinates)-Coordinates(NCoordinates-1);
  
  % average value of right cell
  vm_r=v(1);
  
  % left values, exact
  v_pos=lshap(1:dofs)*v(offset+(1:dofs));
      
  % right values exact
  v_neg=rshap(1:dofs)*v(offset+(1:dofs));
 
  u(offset+(1:dofs))=v(offset+(1:dofs));
  u(offset+2)=minmod(v(offset+2),...
             (vm_m-vm_l)-(v_pos-vm_m-v(offset+2)),...
             (vm_r-vm_m)-(v_neg-vm_m-v(offset+2)),m*h^2);
 
  end % p~=0
  return