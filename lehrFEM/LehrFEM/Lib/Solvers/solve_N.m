function  Un = solve_N(Ln,Sn,gn,NCDofs,varargin)
% SOLVE_N Compute solution in nodal finite element space.
%
%   UN = SOLVE_N(LN,SN,GN,NCDOFS) computes the solution UN in the nodal
%   finite element space.
%
%   Example:
% 
%   Un = solve_N(Ln,Sn,gn,NCDofs)

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland  

  Un = Sn\[Ln; gn];
  Un = Un(NCDofs);
    
return
