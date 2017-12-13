% 2D Wave equation with perfectly matched layer (PML)
%
% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NSTEPS = 2000;                           % Number of time steps
  L = 0.75;                                % Boundary of computational domain
  T = 3;                                   % Final time
  PML = 0.25;                              % Size of PML layer
  SIGMA_0 = 1000;                          % Parameter of absrption profile
  SIGMA_X = @sigma_x_var;                  % Absorption profile in x-direction
  SIGMA_Y = @sigma_y_var;                  % Absorption profile in y-direction 
  RHO = @(x,varargin)ones(size(x,1),1);    % Coefficient function of the wave equation
  GAMMA = @(x,varargin)ones(size(x,1),1);  % Coefficient function of the wave equation
  U0 = @pulse_2D;                          % Initial data
  V0 = @(x,varargin)zeros(size(x,1),1);    % Initial velocity
    
  NFRAMES = 200;                           % Number of frames to be plotted
  FILENAME = 'Wave_2D_PML.avi';            % File name for .avi movie (uncompressed!)
  XLim = (L+PML)*[-1.05 1.05];             % x-axes limits of plot
  YLim = XLim;                             % y-axes limits of plot
  CLim = [-0.2 1.2];                       % Color axes limits of plot
  X_DOM = L*[-1  1  1 -1 -1];              % Boundary of computational domain (x-coordinates)
  Y_DOM = L*[-1 -1  1  1 -1];              % Boundary of computational domain (y-coordinates)
  Z_DOM = CLim(2)*ones(1,5);               % Boundary of computational domain (z-coordinates)
  X_PML = (L+PML)*[-1  1  1 -1 -1];        % Boundary of PML layer (x-coordinates)
  Y_PML = (L+PML)*[-1 -1  1  1 -1];        % Boundary of PML layer (y-coordinates)
  Z_PML = CLim(2)*ones(1,5);               % Boundary of PML layer (z-coordinates)
  
  % Open up FIFO buffer
 
  buf = open(buffer(),'/scratch/users/meury');
  
  % Initialize mesh
  
  Mesh = genMesh_PML(L,PML,0.04);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1; 
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  nElements = size(Mesh.Elements,1);
 
  % Precompute stiffness and mass matrices
  
  dt = T/NSTEPS;
  
  ML = assemMat_LFE(Mesh,@MASS_LFE); 
  [IL,JL,ML_rho] = assemMat_LFE(Mesh, ...
                                @MASS_Weight_LFE, ...
                                P3O3(),RHO);
  [IL,JL,ML_sigma_x] = assemMat_LFE(Mesh, ...
                                    @MASS_Weight_LFE, ...
                                    P3O3(),SIGMA_X,SIGMA_0,L);
  [IL,JL,ML_sigma_y] = assemMat_LFE(Mesh, ...
                                    @MASS_Weight_LFE, ...
                                    P3O3(),SIGMA_Y,SIGMA_0,L);
  
  [IC,JC,MC] = assemMat_P0(Mesh,@MASS_P0);
  [IC,JC,MC_sigma_x] = assemMat_P0(Mesh, ...
                                   @MASS_Weight_P0, ...
                                   P3O3(),SIGMA_X,SIGMA_0,L);
  [IC,JC,MC_sigma_y] = assemMat_P0(Mesh, ...
                                   @MASS_Weight_P0, ...
                                   P3O3(),SIGMA_Y,SIGMA_0,L);
  
  G_x = assemMat_P1P0(Mesh, ...
                      @STIMA_Div_x_P1P0, ...
                      P3O3(),@(x,varargin)ones(size(x,1),1));
  G_x_gamma = assemMat_P1P0(Mesh, ...
                            @STIMA_Div_x_P1P0, ...
                            P3O3(),GAMMA);
  G_y = assemMat_P1P0(Mesh, ...
                      @STIMA_Div_y_P1P0, ...
                      P3O3(),@(x,varargin)ones(size(x,1),1));
  G_y_gamma = assemMat_P1P0(Mesh, ...
                            @STIMA_Div_y_P1P0, ...
                            P3O3(),GAMMA);  

  SL_x = sparse([IL IL],[JL JL],[ML_rho/dt ML_sigma_x/2]);
  SL_y = sparse([IL IL],[JL JL],[ML_rho/dt ML_sigma_y/2]);
  SC_x = sparse([IC IC],[JC JC],[MC/dt MC_sigma_x/2]);
  SC_y = sparse([IC IC],[JC JC],[MC/dt MC_sigma_y/2]);
  
  ML_rho = sparse(IL,JL,ML_rho);
  ML_sigma_x = sparse(IL,JL,ML_sigma_x);
  ML_sigma_y = sparse(IL,JL,ML_sigma_y);
  
  MC = sparse(IC,JC,MC);
  MC_sigma_x = sparse(IC,JC,MC_sigma_x);
  MC_sigma_y = sparse(IC,JC,MC_sigma_y);
  
  G_x_gamma = transpose(G_x_gamma);
  G_y_gamma = transpose(G_y_gamma);
  
  clear('IL','Jl','IC','JC');
  
  % Compute initial data
  
  b = assemLoad_LFE(Mesh,P3O3(),U0);
  b = (ML\b)/2;
  U_x_new = b;
  U_y_new = b;
  
  V_x_new = zeros(nElements,1);
  V_y_new = zeros(nElements,1);
 
  b = assemLoad_LFE(Mesh,P3O3(), ...
                    @(x,varargin)RHO(x,varargin{:}).*V0(x,varargin{:}));
                    
  U = U_x_new + U_y_new;
  buf = push(buf,U);
      
  % Integrate ODE system (dissipative leapfrog)
 
  per = 0;
  progress_bar(per);
  U_x_old = U_x_new;
  U_y_old = U_y_new;
  V_x_old = V_x_new;
  V_y_old = V_y_new;
  for i = 1:NSTEPS
    if(per < floor(100*i/NSTEPS))
      per = floor(100*i/NSTEPS);
      progress_bar(per);
    end  
    
    % Incorporate Dirichlet boundary conditions  
      
    [U_x_new,FreeDofs] = assemDir_LFE(Mesh,-1,@(x,varargin)zeros(size(x,1),1));
    [U_y_new,FreeDofs] = assemDir_LFE(Mesh,-1,@(x,varargin)zeros(size(x,1),1));
    
    % Update x-component of U
            
    rhs = b + (ML_rho*U_x_old)/dt - (ML_sigma_x*U_x_old)/2 - G_x*V_x_old;
    U_x_new(FreeDofs) = SL_x(FreeDofs,FreeDofs)\rhs(FreeDofs);
    
    % Update y-component of U
    
    rhs = b + (ML_rho*U_y_old)/dt - (ML_sigma_y*U_y_old)/2 - G_y*V_y_old;
    U_y_new(FreeDofs) = SL_y(FreeDofs,FreeDofs)\rhs(FreeDofs);
  
    % Save solution to FIFO buffer
    
    U = U_x_new + U_y_new;
    if(rem(i,max(floor(NSTEPS/NFRAMES),1)) == 0)
      buf = push(buf,U);
    end
    
    % Update x-component of V
    
    rhs = (MC*V_x_old)/dt - (MC_sigma_x*V_x_old)/2 + G_x_gamma*U;
    V_x_new = SC_x\rhs;
    
    % Update y-component of V
    
    rhs = (MC*V_y_old)/dt - (MC_sigma_y*V_y_old)/2 + G_y_gamma*U;
    V_y_new = SC_y\rhs;
          
    % Update temporary variables
    
    U_x_old = U_x_new;
    U_y_old = U_y_new;
    V_x_old = V_x_new;
    V_y_old = V_y_new;
    
  end
  
  % Generate movie
  
  if(~isempty(FILENAME))
    Mov = avifile(FILENAME);
  end
  fig = figure('Name','2D Wave equation with PML');
  while(~isempty(buf))
    
    % Extract data from FIFO buffer  
      
    [U,buf] = pop(buf);
    
    % Generate figure
    
    clf;
    patch('Faces', Mesh.Elements, ...
          'Vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) U], ...
          'CData', U, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    hold on;
    plot3(X_DOM,Y_DOM,Z_DOM,'k:', ...
          X_PML,Y_PML,Z_PML,'k-');
    hold off;
    set(gca,'XLim',XLim, ...
            'YLim',YLim, ...
            'CLim',CLim, ...
            'DataAspectRatio',[1 1 1]);
    colorbar;
    drawnow();
    
    % Add frame to movie
    
    if(~isempty(FILENAME))
      F = getframe();
      Mov = addframe(Mov,F);
    end
    
  end
  buf = close(buf);
    
  % Close the movie
  
  if(~isempty(FILENAME))
    Mov = close(Mov);
  end
  
  % Clear previously used memory
  
  clear all;
  close all;
  