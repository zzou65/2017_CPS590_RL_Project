function varargout = plot_P0Cell(U,Mesh)
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

%   Copyright 2006-2006 Patrick Meury & Kah-Ling Sia & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  CHAR_FUN_Handle=@(x,varargin)ones(size(x,1),1);
  
  b_Mesh=refine_BAR(Mesh);
  b_Mesh = add_Edge2Elem(b_Mesh);
  b_Mesh = add_Patches(b_Mesh);
  b_nElements=size(b_Mesh,1);
  b_U=zeros(b_nElements,1);
  nCoordinates=size(Mesh.Coordinates);
  M=assemCochainD_2f(Mesh,CHAR_FUN_Handle,P1O2());
  for i=1:nCoordinates
      b_AdjElements=b_Mesh.AdjElements(i,:);
      b_AdjElements=setdiff(b_AdjElements,0);
      b_n_AdjElements=size(b_AdjElements,2);
      b_U(b_AdjElements)=U(i)/M(i);
  end
  plot_P0(b_U,b_Mesh);

  % Compute axis limits
  
return
