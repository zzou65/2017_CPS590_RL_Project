function err = L2Err_P2f(Mesh,u,QuadRule,FHandle)
% error finite element solution for Primal 2 forms
%
%   
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  CHAR_FUN_Handle=@(x,varargin)ones(size(x,1),1);
  
  M=assemCochain_2f(Mesh,CHAR_FUN_Handle,P1O2());
  u=u./M;
  err=L2Err_PC(Mesh,u,QuadRule,FHandle);
  
return