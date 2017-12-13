% Convergence rates for 1D hpDG discretization schemes.
% Laplace equation

% Copyright 2007-2007 Patrick Meury & Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize mesh
  
  X0 = -1;                    % Left end point of interval
  X1 =  1;                    % Right end point of interval
  NPOINTS = 4;                % Initial set of grid points
  NREFS = 1;                 % Number of grid refinements
  PMAX = 20;                   % Maximal polynomial degree
  PMIN=1;
  UEX = @(x,varargin)4*x.^1+3;  % Interpolated function
  s= 1;                       % switch between: 1 SIP, -1 NIP, 0 IIP  
  alpha=10;                   % penaltyparameter
  % Initialize quadrature rule and shape functions
  
  qr = gauleg(-1,1,ceil(PMAX+1/2));
  Grad_Shap = grad_shap_Leg_1D(qr.x,PMAX);
  Shap = shap_Leg_1D(qr.x,PMAX);
  
  Inn_shap=shap_Leg_1D([-1;1],PMAX);
  Inn_gradshap=grad_shap_Leg_1D([1;-1],PMAX);
  % Compute discretization error
  
  h = zeros(NREFS,1);
  err = zeros(NREFS,PMAX-1);
  for i = 1:NREFS
    for j = PMIN:1:PMAX
  
      % Initialize mesh
  
      Coordinates = X0 + (X1-X0)/(NPOINTS-1)*(0:(NPOINTS-1));
  
      % Initialize polynomial degrees
  
      p = j*ones(1,NPOINTS-1);
      
      % Assemble mass matrices and load vector
  
      %A_vol = assemMat_Vol_hpDG_1D(Coordinates,p,@STIMA_Lapl_hpDG_1D,qr,Shap);
      A_Inn =assemMat_Inn_hpDG_1D(Coordinates,p,@STIMA_Inn_hpDG_1D,Inn_shap,Inn_gradshap,s);
      A_InnPen=assemMat_Inn_hpDG_1D(Coordinates,p,@STIMA_InnPen_hpDG_1D,Inn_shap,Inn_gradshap,alpha);
      A_Bnd=assemMat_Bnd_hpDG_1D(Coordinates,p,Inn_shap,Inn_gradshap,s,alpha);
      M=assemMat_Vol_hpDG_1D(Coordinates,p,@MASS_Vol_hpDG_1D,qr,Shap);
      L=assemLoad_hpDG_1D(Coordinates,p,qr,Shap,UEX);

      %S=A_vol-A_Inn+A_InnPen;
      % Solve the linear system  
      u = (M)\L;
    
      % Compute discretization error
     
      err(i,j+1) = L2Err_hpDG_1D(Coordinates,p,u,qr,Shap,UEX);
      
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
  set(gca,'XScale','log', ...
          'XDir','reverse', ...
          'YScale','log');
  title('{\bf L^2 Error}');
      
  % Clear memory
  
  clear all;