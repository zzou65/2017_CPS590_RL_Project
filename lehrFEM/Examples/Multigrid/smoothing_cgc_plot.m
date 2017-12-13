function [] = smoothing_cgc_plot()
% plot errors for coarse grid correction after various amounts of smoothing
%
%   This code plots the error distribution of a sample problem after 0, 1,
%   2 and 3 steps of Gauss-Seidel smoothing and a subsequent coarse grid
%   corection.  The coarse grid correction seems to have little effect
%   without the presmoothing; however, there is a clear decrease in the
%   error when the coarse grid correction is applied after smoothing.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize constants
 
  ref = [3 4];                                   % refinements of initial mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  its = 3;                                       % number of smoothing iterations
  smoother = @gs_smooth;                         % multigrid smoother
  
  % generate multigrid data structure
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',Mesh,'ref',ref);
  mg_data = mg_stima(mg_data,'f',f_handle,'gd',gd_handle);
  mg_data = mg_error(mg_data,'energy',false,'iter',false);
  
  % plot meshes
  
  plot_Mesh(Mesh,'as');
  plot_Mesh(mg_data{1}.mesh,'as');
  plot_Mesh(mg_data{2}.mesh,'as');
  
  % calculate exact solution
  
  U_bd = mg_data{2}.u_bd;
  U_ex = U_bd;
  U_ex(mg_data{2}.dofs) = mg_data{2}.A\mg_data{2}.b;
  
  % define oscillatory initial guess
  
  u0 = @(x) cos(0.5*pi*x(:,1)).*cos(0.5*pi*x(:,2))+sin(42*pi*x(:,1)).*sin(57*pi*x(:,2))+x(:,1).^2+x(:,2).^2;
  U = repmat(U_bd,1,its+1);
  U(mg_data{2}.dofs,1) = u0(mg_data{2}.mesh.Coordinates(mg_data{2}.dofs,:));
  
  
  % initialize errors
  
  Err = zeros(size(U));
  
  % carry out multigrid v-cycles with various amounts of (pre)smoothing
  
  for i=1:its+1
    mg_data = mg_smooth(mg_data,'m',[i-1,0],'smoother',smoother);
    U(mg_data{2}.dofs,i) = mg(U(mg_data{2}.dofs,1),mg_data,0,1);
    Err(:,i) = U(:,i)-U_ex;
  end
  
  % plot errors
  
  for i=1:its+1
    plot_LFE(Err(:,i),mg_data{2}.mesh);
%     set(gca,'CameraPosition',[1,2,1]);
%     grid on;
%     axis([-1,1,-1,1,-1,1])
    set(gca,'CameraPosition',[1,2,1.5]);
    grid on;
    axis([-1,1,-1,1,-0.5,1])
    xlabel('\bf x');
    ylabel('\bf y');
    zlabel('\bf error');
    title(sprintf('%s Error After %.0f Smoothing Step(s) and CGC','\bf',i-1));
  end
  
return