% Convergence rate for normal and smart flux computation using linear and
% quadratic finite elements for the heat equation with Dirichlet boundary
% conditions on the unit square. This script generates the following .eps
% files:
%  Pipe_flux_meshwidth,  Pipe_smartflux_meshwidth.eps.

% Copyright 2005-2005 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% 8092-Zurich, Switzerland

  % Initialize constants
  
  R_OUT = 1;                                      % Radius of outer circle
  R_IN = 0.5;                                     % Radius of inner circle
  C_IN = 60;                                      % Temperature inside the tube
  C_OUT = 10;                                     % Temperature outside the tube  
  SIGMA = 1;                                      % Heat conductivity                                                                                                      
  H0 = 0.25;                                      % Initial mesh width
  NREFS = 5;                                      % Number of red refinements
  NGAUSS = 3;                                     % Number of Gauss points
  P_HANDLE = @(x,varargin)-2/(R_IN^2-R_OUT^2)*x;  % Gradient of cut-off function   
  
  % Initialize coefficient functions
  
  C1 = (C_OUT-C_IN)/(log(R_OUT)-log(R_IN));
  C2 = (log(R_OUT)*C_IN-log(R_IN)*C_OUT)/(log(R_OUT)-log(R_IN));
  S_HANDLE = @(x,varargin)SIGMA*ones(size(x,1),1);
  F_HANDLE = @(x,varargin)zeros(size(x,1),1);
  GD_HANDLE = @(x,varargin)C1*log(sqrt(x(:,1).^2+x(:,2).^2))+C2*ones(size(x,1),1);
  JEX = 2*pi*SIGMA*C1;
  
  % Initialize mesh
  
  DHandle = inline(['dist_diff(dist_circ(x,[0 0],' num2str(R_OUT) ...
                    '),dist_circ(x,[0 0],' num2str(R_IN) '))'],'x');
  Mesh = init_Mesh([-R_OUT -R_OUT; R_OUT R_OUT],H0,DHandle,@h_uniform,[],0);
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);  
  MidPts = 1/2*(Mesh.Coordinates(Mesh.Edges(Loc,1),:) + ...
                Mesh.Coordinates(Mesh.Edges(Loc,2),:));
  Loc_in = Loc(sqrt(sum((MidPts).^2,2)) > (R_OUT+R_IN)/2);
  Loc_out = Loc(sqrt(sum((MidPts).^2,2)) < (R_OUT+R_IN)/2);
  Mesh.BdFlags(Loc_in) = -1;
  Mesh.BdFlags(Loc_out) = -2;
  
  % Compute discretization error on a series of meshes
    
  h = zeros(1,NREFS);
  JErr_LFE = zeros(1,NREFS);
  JErr_QFE = zeros(1,NREFS);
  JErr_smart_LFE = zeros(1,NREFS);
  JErr_smart_QFE = zeros(1,NREFS);
  for i = 1:NREFS
    
    % Refine the mesh 
      
    Mesh = refine_REG(Mesh,DHandle);
    Mesh = add_Edge2Elem(Mesh);
   
    % Assemble stiffness matrix and load vector
   
    A_LFE = assemMat_LFE(Mesh,@STIMA_Heat_LFE,P7O6(),S_HANDLE);
    L_LFE = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
        
    A_QFE = assemMat_QFE(Mesh,@STIMA_Heat_QFE,P7O6(),S_HANDLE);
    L_QFE = assemLoad_QFE(Mesh,P7O6(),F_HANDLE);
    
    % Incorporate Dirichlet boundary conditions
    
    [U_LFE,FreeDofs_LFE] = assemDir_LFE(Mesh,-1:-1:-2,GD_HANDLE);
    L_LFE = L_LFE - A_LFE*U_LFE;
    
    [U_QFE,FreeDofs_QFE] = assemDir_QFE(Mesh,-1:-1:-2,GD_HANDLE);
    L_QFE = L_QFE - A_QFE*U_QFE;
    
    % Solve the linear system
  
    U_LFE(FreeDofs_LFE) = A_LFE(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    U_QFE(FreeDofs_QFE) = A_QFE(FreeDofs_QFE,FreeDofs_QFE)\L_QFE(FreeDofs_QFE);
        
    % Compute fluxes
    
    JErr_LFE(i) = abs(Flux_LFE(U_LFE,-1,Mesh,gauleg(0,1,NGAUSS),S_HANDLE)-JEX);
    JErr_QFE(i) = abs(Flux_QFE(U_QFE,-1,Mesh,gauleg(0,1,NGAUSS),S_HANDLE)-JEX);
    JErr_smart_LFE(i) = abs(SmartFlux_LFE(U_LFE,Mesh,P_HANDLE,P7O6(),S_HANDLE)-JEX);
    JErr_smart_QFE(i) = abs(SmartFlux_QFE(U_QFE,Mesh,P_HANDLE,P7O6(),S_HANDLE)-JEX);
    
    h(i) = get_MeshWidth(Mesh);
    
  end
  
  save('PipeFlowErrs.mat','h','JErr_LFE','JErr_QFE','JErr_smart_LFE','JErr_smart_QFE');
  
  % Plot convergence rate of boundary flux for standard flux and linear FEM 
  fig = figure;
  plot(h,JErr_LFE,'b-',h,JErr_LFE,'k+'); grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  ax = axis; axis([0.01 0.5 ax(3) ax(4)]); 
  title('{\bf Convergence: standard boundary flux, linear FE}');
  xlabel('{\bf Mesh width [log scale]}');
  ylabel('{\bf Boundary flux [log scale]}');
  p = polyfit(log(h),log(JErr_LFE),1);
  add_Slope(gca,'South',p(1));
  print('-depsc2','Bdflux_meshwidth.eps');
  
  % Plot error in both fluxes 
  fig = figure;
  plot(h,JErr_LFE,'b-',h,JErr_smart_LFE,'r-',h,JErr_LFE,'k+',h,JErr_smart_LFE,'k*'); 
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  ax = axis; axis([0.01 0.5 ax(3) ax(4)]);
  legend('Unstable flux J','Stable flux J^{\ast}','Location','East');
  title('{\bf Convergence: stable boundary flux, linear FE}');
  xlabel('{\bf Mesh width [log scale]}');
  ylabel('{\bf Boundary flux [log scale]}');
  p = polyfit(log(h),log(JErr_LFE),1);
  add_Slope(gca,'NorthWest',p(1));
  p = polyfit(log(h),log(JErr_smart_LFE),1);
  add_Slope(gca,'South',p(1));
  print('-depsc2','StabBdflux_meshwidth.eps');
  
  % Plot out normal flux convergence rate and add slope triangles to the plot
  
  fig = figure;
  plot(h,JErr_LFE,'r-', ...
       h,JErr_QFE,'b-', ...
       h,JErr_LFE,'k+', ...
       h,JErr_QFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  ax = axis; axis([0.01 0.5 ax(3) ax(4)]);
  title('{\bf Convergence: instable boundary flux}');
  xlabel('{\bf Mesh width [log scale]}');
  ylabel('{\bf Flux [log scale]}');
 
  legend('Linear FE','Quadratic FE','Location','East');
  p = polyfit(log(h),log(JErr_LFE),1);
  add_Slope(gca,'South',p(1));
  p = polyfit(log(h),log(JErr_QFE),1);
  add_Slope(gca,'NorthWest',p(1));
  
  print('-depsc','Pipe_flux_meshwidth.eps');
  % !gv Pipe_flux_meshwidth.eps &
  
  % Plot out smart flux convergence rates and add slope triangles to the plot
  
  fig = figure;
  plot(h,JErr_smart_LFE,'r-', ...
       h,JErr_smart_QFE,'b-', ...
       h,JErr_smart_LFE,'k+', ...
       h,JErr_smart_QFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  ax = axis; axis([0.01 0.5 ax(3) ax(4)]);
  title('{\bf Convergemce: stable boundary flux}');
  xlabel('{\bf Mesh width [log scale]}');
  ylabel('{\bf Flux [log scale]}');
   
  legend('Linear FE','Quadratic FE','Location','East');
  p = polyfit(log(h),log(JErr_smart_LFE),1);
  add_Slope(gca,'South',p(1));
  p = polyfit(log(h),log(JErr_smart_QFE),1);
  add_Slope(gca,'NorthWest',p(1));
   
  print('-depsc','Pipe_sflux_meshwidth.eps');
  % !gv Pipe_sflux_meshwidth.eps &
  
  % Clear memory
  
  clear all;
  