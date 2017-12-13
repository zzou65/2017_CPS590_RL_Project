function Lhu = assemRKDG_update_Inn(Coordinates,u,p,Inn_shap,numflowhandle,varargin)
% RKDG_update_Inn
%
% update for cell boundary terms for RKDG
%
%
%   Coordinates: n+1 Meshpoints
%             u: solution in last time step
%             p: polynomial approxiamtion each of the n cells
%         lshap: basisfunctions evaluated left boundary
%         rshap: basisfunctions evaluated right boundary
%    flowhandle: flux, numericalflux 
%      Quadrule: 1D quadrature rule 
%   Example:
%   
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NCoordinates=size(Coordinates,2);
  Dofs=sum(p+1);
  lshap=Inn_shap(2,:);
  rshap=Inn_shap(1,:);
  % Alocate memory
  Lhu=zeros(Dofs,1);
  
  % left and right boundary of domain
  % left value
  offset=sum(p(1:(NCoordinates-2))+1); 
  dofs_l=p(NCoordinates-1)+1;
  u_l=lshap(1:dofs_l)*u(offset+(1:dofs_l));
      
  % right value
  dofs_r=p(1)+1;
  u_r=rshap(1:dofs_r)*u(1:dofs_r);
      
  % calculate limiter here ??????????
      
  % calculate numerical flux
  h=numflowhandle(u_l,u_r);
      
  Lhu(offset+(1:dofs_l))=h*[lshap(1:dofs_l)];
  Lhu(1:dofs_r)=h*[-rshap(1:dofs_r)];
      
  offset=0;                    
  % loop over boundary of cells
  for i=2:(NCoordinates-1)
      
      % left value
      dofs_l=p(i-1)+1;
      u_l=lshap(1:dofs_l)*u(offset+(1:dofs_l));
      
      % right value
      dofs_r=p(i)+1;
      u_r=rshap(1:dofs_r)*u(offset+dofs_l+(1:dofs_r));
      
      % calculate limiter here ????
      
      % calculate numerical flux
      h=numflowhandle(u_l,u_r);
      
      Lhu(offset+(1:dofs_l+dofs_r))=Lhu(offset+(1:dofs_l+dofs_r))+...
          h*[lshap(1:dofs_l)' ; -rshap(1:dofs_r)'];
      
      offset=offset+dofs_l;
  end
  
  return