% 1D Wave equation with perfectly matched layer (PML)
%
% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  L = 2;                                                                                                % Length of the interval
  T = 1.5;                                                                                             % Final time
  C = 2;                                                                                                % Speed of sound
  ALPHA = 0.75;                                                                                  % Right travelling wave
  BETA = 0.25;                                                                                     % Left travelling wave
  U0 = @(x,varargin)(ALPHA+BETA)*pulse_1D(x,varargin{:});        % Initial data
  V0 = @(x,varargin)(ALPHA-BETA)*C*dpulse_1D(x,varargin{:});  % Initial velocity
  MU = @(x,varargin)C^2*ones(size(x,1),1);                                    % Coefficient of the wave equation
  RHO = @(x,varargin)ones(size(x,1),1);                                          % Coefficient of the wave equation  
  NPTS = 1000;                                                                                    % Number of points
  NSTEPS = 10000;                                                                              % Number of time steps

  PML = 0.2;                                                                                         % Length of PML layer
  SIGMA_0 = 100;                                                                                % Scaling parameter for absorption profile
  SIGMA = @sigma_const_1D;                                                            % Absorption profile
    
  NFRAMES = 200;                                                                             % Number of frames
  FILENAME = 'Wave_1D_PML.avi';                                                    % Filename of .avi file
  XLim = [-(L+PML) L+PML];                                                             % X-axes limits
  YLim = [-0.05 1.05];                                                                      % Y-axes limits
  
  % Initialize mesh
 
  dx = (XLim(2)-XLim(1))/(NPTS-1);
  Coordinates = transpose(XLim(1):dx:XLim(2));
  
  % Precompute matrices
  
  dt = T/NSTEPS; 
  QuadRule = gauleg(0,1,4);
  
  MC = assemMat_P0_1D(Coordinates,@MASS_P0_1D);
  ML = assemMat_P1_1D(Coordinates,@MASS_P1_1D);
  ML_sigma = assemMat_P1_1D(Coordinates,@MASS_Weight_P1_1D, ...
                            QuadRule,SIGMA,SIGMA_0,L);
  MC_sigma = assemMat_P0_1D(Coordinates,@MASS_Weight_P0_1D, ...
                            QuadRule,SIGMA,SIGMA_0,L);  
  ML_rho = assemMat_P1_1D(Coordinates,@MASS_Weight_P1_1D, ...
                          QuadRule,RHO);
  G = assemMat_P1P0_1D(Coordinates,@STIMA_Div_P1P0_1D, ...
                       QuadRule,@(x,varargin)1);
  G_mu = assemMat_P1P0_1D(Coordinates,@STIMA_Div_P1P0_1D, ...
                          QuadRule,MU);
  
  G_mu = transpose(G_mu);
  S1 = ML_rho/dt + ML_sigma/2;
  S2 = ML_rho/dt - ML_sigma/2;
  S3 = MC/dt + MC_sigma/2;
  S4 = MC/dt - MC_sigma/2;
  
  % Open up FIFO buffer
 
  buf = open(buffer(),'/home/hheumann');
  
  % Compute initial data
 
  U_old = ML\assemLoad_P1_1D(Coordinates,QuadRule,U0);
  V_old = zeros(NPTS-1,1);
  buf = push(buf,U_old);
  
  b = assemLoad_P1_1D(Coordinates,QuadRule, ...
                      @(x,varargin)RHO(x,varargin{:}).*V0(x,varargin{:}));
  
  % Integrate ODE system (dissipative leapfrog scheme, CFL condition)
  
  FreeDofs = 2:(NPTS-1);
  for i = 1:NSTEPS
    
    % Compute right hand side load data
    
    rhs = b + S2*U_old - G*V_old;
    
    % Incorporate zero Dirichlet boundary data
    
    U_new = zeros(NPTS,1);
    
    % Compute new value for U
    
    U_new(FreeDofs) = S1(FreeDofs,FreeDofs)\rhs(FreeDofs);
    if(rem(i,max(floor(NSTEPS/NFRAMES),1)) == 0)
      buf = push(buf,U_new);
    end
    
    % Compute new value for V
    
    rhs = S4*V_old + G_mu*U_new;
    V_new = S3\rhs;

    % Update old values for U and V
    
    U_old = U_new;
    V_old = V_new;
    
  end
  
  % Generate movie
  
  if(~isempty(FILENAME))
    Mov = avifile(FILENAME);
  end
  fig = figure('Name','1D Wave equation with PML');
  while(~isempty(buf))
    
    % Extract data from FIFO buffer  
      
    [U,buf] = pop(buf);
    
    % Generate figure
    
    clf;
    plot(Coordinates,U,'r-', ...
         [-L -L],YLim,'k--', ...
         [ L  L],YLim,'k--');
    set(gca, ...
        'XLim',XLim, ...
        'YLim',YLim);
    
    % Add frame to movie
    
    if(~isempty(FILENAME))
      F = getframe();
      Mov = addframe(Mov,F);
    end
    
  end
    
  % Close the movie
  
  if(~isempty(FILENAME))
    Mov = close(Mov);
  end
  
  % Clear memory and closes figure
  
  clear all;
  close all;
  