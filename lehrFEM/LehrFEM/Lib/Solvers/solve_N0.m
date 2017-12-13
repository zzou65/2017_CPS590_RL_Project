function Un = solve_N0(Ln,Sn,varargin)
% SOLVE_N0 Compute solution in nodal finite element space.
%
%   UN = SOLVE_N0(LN,SN) computes the solution N0 in the nodal finite
%   element space with zero tangential component on the boundary.
%
%   Example:
% 
%   Un = solve_N0(Ln,Sn)

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland  

  Un = Sn\Ln;
    
return
    
