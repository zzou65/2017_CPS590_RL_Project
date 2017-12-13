function u = SlopeLim_lr(Coordinates,v,p,Inn_shap, shap, QuadRule,m,varargin)
% less restictive slope limiter according to Cockburn & Shu
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
%   v_cell=shap(:,1:p(NElements)+1)*v(sum(p(1:NElements-1)+1)+(1:p(NElements)+1));
%   vm_l = sum(QuadRule.w.*v_cell)/2; %*h/h
  
  % middle cell average
  vm_m=v(1);
%   v_cell=shap(:,1:p(1)+1)*v(1:p(1)+1);
%   vm_m = sum(QuadRule.w.*v_cell)/2;
  
  offset=0;                    
  % loop over cells
  for i=1:(NElements-1)
     
      % Cellsize and dofs of middle cell
      dofs=p(i)+1;
      h=Coordinates(i+1)-Coordinates(i); 
      
      % right cell average
      vm_r=v(offset+dofs+1);
%       v_cell=shap(:,1:p(i+1)+1)*v(offset+dofs+(1:p(i+1)+1));
%       vm_r = sum(QuadRule.w.*v_cell)/2;
      
      % limiter calculation
      
      % left values, exact and limited
      v_pos=lshap(1:dofs)*v(offset+(1:dofs));
      u_pos=vm_m-minmod(vm_m-v_pos,vm_m-vm_l,vm_r-vm_m,m*h^2);
      
      % right values exact and limited
      v_neg=rshap(1:dofs)*v(offset+(1:dofs));
      u_neg=vm_m+minmod(v_neg-vm_m,vm_m-vm_l,vm_r-vm_m,m*h^2);
      
      if ((u_neg==v_neg) && (u_pos==v_pos))
          u(offset+(1:dofs))=v(offset+(1:dofs));
      else
         u(offset+(1:dofs))=0;
         u(offset+1)=v(offset+1);
         u(offset+2)=minmod(v(offset+2),(vm_m-vm_l),(vm_r-vm_m),m*h^2);
      end
      offset=offset+dofs;
      vm_l=vm_m;
      vm_m=vm_r;
      
  end
  dofs=p(NElements)+1;
  h=Coordinates(NCoordinates)-Coordinates(NCoordinates-1);
  
  % average value of right cell
  vm_r=v(1);
%   v_cell=shap(:,1:p(1)+1)*v(1:p(1)+1);
%   vm_r = sum(QuadRule.w.*v_cell)/2;
  
  % left values, exact and limited
  v_pos=lshap(1:dofs)*v(offset+(1:dofs));
  u_pos=vm_m-minmod(vm_m-v_pos,vm_m-vm_l,vm_r-vm_m,m*h^2);
      
  % right values exact and limited
  v_neg=rshap(1:dofs)*v(offset+(1:dofs));
  u_neg=vm_m+minmod(v_neg-vm_m,vm_m-vm_l,vm_r-vm_m,m*h^2);
      
  if (u_neg==v_neg) && (u_pos==v_pos)
      u(offset+(1:dofs))=v(offset+(1:dofs));
  else
      u(offset:(1+dofs))=0;
      u(offset+1)=v(offset+1);
      u(offset+2)=minmod(v(offset+2),(vm_m-vm_l),(vm_r-vm_m),m*h^2);
  end
  end % p~=0
  return