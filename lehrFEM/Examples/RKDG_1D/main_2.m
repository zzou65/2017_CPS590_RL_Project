% Convergence rates for 1D hpDG discretization schemes.

% Copyright 2007-2007 Patrick Meury & Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize mesh
  
  X0 = -1;                    % Left end point of interval
  X1 =  1;                    % Right end point of interval
  NPOINTS = 2;                % Initial set of grid points
  NREFS = 15;                 % Number of grid refinements
  PMAX = 5;                   % Maximal polynomial degree
  UEX = @(x,varargin)exp(x);  % Interpolated function
    
  % Initialize quadrature rule and shape functions
  
  qr = gauleg(-1,1,ceil(PMAX+1+1/2));
  Shap = shap_Leg_1D(qr.x,PMAX);
  
  % Compute discretization error
  
  h = zeros(NREFS,1);
  err = zeros(NREFS,PMAX+1);
  for i = 1:NREFS
    for j = 0:1:PMAX
  
      % Initialize mesh
  
      Coordinates = X0 + (X1-X0)/(NPOINTS-1)*(0:(NPOINTS-1));
  
      % Initialize polynomial degrees
  
      p = j*ones(1,NPOINTS-1);
      
      % Assemble mass matrices and load vector
  
      M = assemMat_Vol_hpDG_1D(Coordinates,p,@MASS_Vol_hpDG_1D,qr,Shap);
      L = assemLoad_hpDG_1D(Coordinates,p,qr,Shap,UEX);

      % Solve the linear system
  
      u = M\L;
     % plot_hpDG_1D(Coordinates,p,u,@shap_Leg_1D);
      % Compute discretization error
     
      err(i,j+1) = L2Err_hpDG_1D(Coordinates,p,u,qr,Shap,UEX)
     
    end

    % Refine grid
    
    h(i) = (X1-X0)/(NPOINTS-1);
    NPOINTS = 2*NPOINTS;
    
  end  
  
  % Generate figure
  
  fig = figure('Name','L2 Error');  
  hold on;
  for i = 0:1:PMAX
    plot(h,err(:,i+1),'r-x');
  end
  hold off;
  set(gca, ...
          'XDir','reverse', ...
          'YScale','log','XScale','log');
  title('{\bf L^2 Error}');
      
  % Clear memory
  
  clear all;
  