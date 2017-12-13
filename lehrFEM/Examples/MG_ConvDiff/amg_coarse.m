function [] = amg_coarse()
% plot 2-level amg coarse grid points for convection-diffusion problems

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%% linear flow

  % define parameters

  refs = 4;
  phi = pi/8;
  v_handle = @(x,varargin) ones(size(x,1),1)*[sin(phi),cos(phi)];
  c = 1e-3;
  
  % initialize and refine mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh = add_Edges(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  
  for k=1:refs
    Mesh = refine_REG(Mesh);
  end
  
  Mesh = rmfield(Mesh,'BdFlags');
  [inflow,outflow,neutral] = get_inflow_outflow(Mesh,v_handle);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(inflow) = -1;
  Mesh.BdFlags(outflow) = -2;
  Mesh.BdFlags(neutral) = -2;
  
  % determine free degrees of freedom
  
  Loc = find(Mesh.BdFlags==-1);
  DDofs = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
  FDofs = setdiff(1:size(Mesh.Coordinates,1),DDofs);
  
  % compute stiffness matrix
  
  A = c*assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  M = assemMat_MassZeroD(Mesh);
  D_fd = assemMat_LFE(Mesh,@LIEUP_FD,v_handle);

  A_fd = A+M*D_fd;
  
  A_fd = A_fd(FDofs,FDofs);
  
  % construct amg data structure
  
  AMGOpt = AMGDefaultOptions;
  AMGOpt.levsmax = 2;

  AMGData = AMGSetup(A_fd,AMGOpt);
  
  % plot coarse grid points
  
  [fig,h] = plot_AMG_coarse(Mesh,AMGData,FDofs);
  
  % plot velocity field
  
  midpt = 0.5*(Mesh.Coordinates(Mesh.Edges(:,1),:)+Mesh.Coordinates(Mesh.Edges(:,2),:));
  v = v_handle(midpt);
  h1 = quiver(midpt(:,1),midpt(:,2),v(:,1),v(:,2),'k');
  
  h = [h,h1];
  legend(h,'coarse nodes','flow','Location','NorthEastOutside');
  
  
%% circular flow

  % define parameters

  refs = 4;
  v_handle = @vel_circ;
  c = 1e-3;
  
  % initialize and refine mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh = add_Edges(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  
  for k=1:refs
    Mesh = refine_REG(Mesh);
  end
  
  Mesh = rmfield(Mesh,'BdFlags');
  [inflow,outflow,neutral] = get_inflow_outflow(Mesh,v_handle);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(inflow) = -1;
  Mesh.BdFlags(outflow) = -2;
  Mesh.BdFlags(neutral) = -2;
  
  % determine free degrees of freedom
  
  Loc = find(Mesh.BdFlags==-1);
  DDofs = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
  FDofs = setdiff(1:size(Mesh.Coordinates,1),DDofs);
  
  % compute stiffness matrix
  
  A = c*assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  M = assemMat_MassZeroD(Mesh);
  D_fd = assemMat_LFE(Mesh,@LIEUP_FD,v_handle);

  A_fd = A+M*D_fd;
  
  A_fd = A_fd(FDofs,FDofs);
  
  % construct amg data structure
  
  AMGOpt = AMGDefaultOptions;
  AMGOpt.levsmax = 2;

  AMGData = AMGSetup(A_fd,AMGOpt);
  
  % plot coarse grid points
  
  [fig,h] = plot_AMG_coarse(Mesh,AMGData,FDofs);
  
  % plot velocity field
  
  midpt = 0.5*(Mesh.Coordinates(Mesh.Edges(:,1),:)+Mesh.Coordinates(Mesh.Edges(:,2),:));
  v = v_handle(midpt);
  h1 = quiver(midpt(:,1),midpt(:,2),v(:,1),v(:,2),'k');
  
  h = [h,h1];
  legend(h,'coarse nodes','flow','Location','NorthEastOutside');
  