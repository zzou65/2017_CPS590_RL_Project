function varargout = plot_P0Simplex(U,Mesh)
% PLOT_P0 Plot finite element solution.
%
%   PLOT_P0(U,MESH) generates a plot of the finite element solution U on
%   the mesh MESH.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the
%                mesh.
%
%   H = PLOT_P0(U,MESH) also returns the handle to the figure.
%
%   Example:
%
%   plot_P0(U,MESH);

%   Copyright 2006-2006 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  CHAR_FUN_Handle=@(x,varargin)ones(size(x,1),1);
  
  M=assemCochain_2f(Mesh,CHAR_FUN_Handle,P1O2());
  U=U./M;
  plot_P0(U,Mesh)
  
return
