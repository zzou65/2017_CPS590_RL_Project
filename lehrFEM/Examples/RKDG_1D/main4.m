% Convergence rates for 1D hpDG discretization schemes.
% Laplace equation

% Copyright 2007-2007 Patrick Meury & Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize mesh
  
  X0 =-0.5;                                       % Left end point of interval
  X1 = 1;                                       % Right end point of interval
  NPOINTS = 3;                                 % Initial set of grid points
  NREFS =5 ;                                    % Number of grid refinements
  PMAX =5   ;                                     % Maximal polynomial degree
  PMIN=1;
  UEX = @(x,varargin)sin(pi*x);                  % Interpolated function
  F = @(x,varargin)sin(pi*x)*pi^2+cos(pi*x)*pi;  % Right hand-side source term
  GD = @(x,varargin)sin(pi*x);                   % Dirichlet boundary data
  Vhandle = @(x,varargin)ones(size(x,2));                   % Dirichlet boundary data
  s= 1;                                         % switch between: 1 SIP, -1 NIP, 0 IIP  
  alpha=10;                                      % penaltyparameter
  % Initialize quadrature rule and shape functions
  
  qr = gauleg(-1,1,ceil(PMAX+1/2));
  error_qr=gauleg(-1,1,100);
  Grad_Shap = grad_shap_Leg_1D(qr.x,PMAX);
  Shap = shap_Leg_1D(qr.x,PMAX);
  error_Shap = shap_Leg_1D(error_qr.x,PMAX);
  
  Inn_shap=shap_Leg_1D([-1;1],PMAX);
  Inn_gradshap=grad_shap_Leg_1D([-1;1],PMAX);
%   
%   Grad_Shap = grad_shap_DGLFE_1D(qr.x);
%   Shap = shap_DGLFE_1D(qr.x);
%   error_Shap = shap_DGLFE_1D(error_qr.x);
%   
%   Inn_shap=shap_DGLFE_1D([-1;1]);
%   Inn_gradshap=grad_shap_DGLFE_1D([-1;1]);
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
  
      A_vol = assemMat_Vol_hpDG_1D(Coordinates,p,@STIMA_Vol_Lapl_hpDG_1D,qr,Grad_Shap);
      C_vol = assemMat_Vol_hpDG_1D(Coordinates,p,@STIMA_Vol_Conv_hpDG_1D,qr,Shap,Grad_Shap,Vhandle);
      A_Inn =assemMat_Inn_hpDG_1D(Coordinates,p,@STIMA_Inn_Lapl_hpDG_1D,Inn_shap,Inn_gradshap,s);
      C_Inn =assemMat_Inn_hpDG_1D(Coordinates,p,@STIMA_Inn_Conv_hpDG_1D,Inn_shap,Vhandle);
      A_InnPen=assemMat_Inn_hpDG_1D(Coordinates,p,@STIMA_InnPen_hpDG_1D,Inn_shap,Inn_gradshap,alpha);
      A_Bnd=assemMat_Bnd_hpDG_1D(Coordinates,p,Inn_shap,Inn_gradshap,s,alpha,Vhandle);
      M=assemMat_Vol_hpDG_1D(Coordinates,p,@MASS_Vol_hpDG_1D,qr,Shap);
      G=assemDir_hpDG_1D(Coordinates,p,GD,Inn_shap,Inn_gradshap,s,alpha,Vhandle); 
      L=assemLoad_hpDG_1D(Coordinates,p,qr,Shap,F);
      U=assemLoad_hpDG_1D(Coordinates,p,qr,Shap,UEX);
      
      L=L-G;
      S=A_vol-C_vol-A_Inn+C_Inn+A_InnPen+A_Bnd;
  
      % Solve the linear system  
      u = S\L;
      
%       plot_hpDG_1D(Coordinates,p,u,@shap_DGLFE_1D);
%       plot_hpDG_1D(Coordinates,p,u,@shap_DGLFE_1D);
        plot_hpDG_1D(Coordinates,p,u,@shap_Leg_1D);
%        plot_hpDG_1D(Coordinates,p,abs(u-M\U),@shap_Leg_1D);
      % Compute discretization error
     
      err(i,j+1) = L2Err_hpDG_1D(Coordinates,p,u,error_qr,error_Shap,UEX);
      
    end

    % Refine grid
    
    h(i) = (X1-X0)/(NPOINTS-1);
    NPOINTS = 2*NPOINTS;
    
  end  
  
  % Generate figure
  
  fig = figure('Name','L2 Error');  
  hold on;
 % plot(PMIN:PMAX,err(PMIN:PMAX),'r-x')
  for i = 0:1:PMAX
    plot(h,err(:,i+1),'r-x');
  end
  hold off;
  set(gca,'XScale','log', ...
          'XDir','reverse', ...
          'YScale','log');
  title('{\bf L^2 Error}');
      
  % Clear memory
  
  %clear all;