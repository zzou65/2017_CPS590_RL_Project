function Mloc = STIMA_GradContr(Vertices,ElemInfo,V_HANDLE,QuadRule,varargin)
% STIMA_GradContr conforming discretisation leads to zero stiffnessmatrix;
%
%   MLOC = STIMA_GradContr(VERTICES,...) computes ... matrix using 
%   Whitney 1-forms finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
% 
%   ElemInfo (not used)
%
%   V_HANDLE 
%   Example:
%

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

Mloc=zeros(3,3);