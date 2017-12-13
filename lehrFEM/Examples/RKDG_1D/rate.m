% Run script for linear advection problem.

% Copyright 2006-2007 Patrick Meury Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland
  
  % Initialize constants
  X0 =0;                 % Left end point of interval
  X1 =2;                 % Right end point of interval
  NPOINTS = 10.^[1:4];           % Initial set of grid points
  P=[0:3];                    % ploynomial degree 
  CFL=0.9*1./(2*P+1);
% CFL=0.1;
  T =20;                   % Final time
 
  U0 = @(x,t)1+0.2*sin(pi*x);                 % Initial data
  Grad_U0 = @(x,t)0.2*pi*cos(pi*x);                 % Initial data
  U0 = @(x,t)(x-2).*(x);                 % Initial data
  Grad_U0 = @(x,t)x-2+x;                 % Initial data
  G = @(x,t)zeros(size(x,1),1);                 % Right hand side load data
  FLUX = @(x) x;                                % Flux function
  NUMFLUX = @(x,y) x;                           % Numerical flux function

  % container for error calculation
  err=zeros(size(NPOINTS,2),size(P,2));
  h=err;
  
  % quadrature rule and shape functions at quadrature point and boundary
  for q=1:size(P,2)
      
      % quadrature rule and shape functions at quadrature point and boundary
      QuadRule=gauleg(-1,1,ceil(P(q)+3/2));
      Shap=shap_Leg_1D(QuadRule.x,P(q));
      Grad_Shap=grad_shap_Leg_1D(QuadRule.x,P(q));
      Inn_shap=shap_Leg_1D([-1;1],P(q));
      Inn_gradshap=grad_shap_Leg_1D([-1;1],P(q));
  
      % error calculation
      qr=gauleg(-1,1,P(q)+50);
      err_Shap=shap_Leg_1D(qr.x,P(q));
  
      for s=1:size(NPOINTS,2);
      
          % Initialize mesh
      
          Coordinates = X0 + (X1-X0)/(NPOINTS(s)-1)*(0:(NPOINTS(s)-1));
          dx=min(Coordinates(2:end)-Coordinates(1:end-1));

          % Initialize polynomial degrees

          p=P(q)*ones(1,NPOINTS(s)-1);
  
          %Assemble mass matrix and extract diagonal
          M=assemMat_Vol_hpDG_1D(Coordinates,p,@MASS_Vol_hpDG_1D,QuadRule,Shap);
          M=spdiags(M,0);
  
          % Compute right data
  
          L = assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,U0,0);
          U = L./M;
          
          L = assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,U0,0);
          U = L./M;
          
          L_g = assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,Grad_U0,0);
     
  
          % Integrate ODE system
            
          Bvol=assemRKDG_update_Vol(Coordinates,U,p,Shap,Grad_Shap,FLUX,QuadRule);
          Binn=assemRKDG_update_Inn(Coordinates,U,p,Inn_shap,NUMFLUX);
          L_p = Bvol-Binn;
          
%            plot_hpDG_1D(Coordinates,p,L_p./M,@shap_Leg_1D)
%            plot_hpDG_1D(Coordinates,p,-L_g./M,@shap_Leg_1D)
             
          % calculate error
    
          err(s,q)=L2Err_hpDG_1D(Coordinates,p,-L_p./M,qr,err_Shap,Grad_U0)
          h(s,q)=dx;
      end % end of meshsize iteration
  end % end of polynomial degree iteration
  
  
  fig = figure('Name','Discretization error');
  hold on;
  ps=size(NPOINTS,2);
  ph=zeros(size(P,2),2);
  for i=1:size(P,2)
    plot(h(:,i),err(:,i),'r+-');
    ph(i,:) = polyfit(log(h(ps-3:ps,i)),log(err(ps-3:ps,i)),1);
  end
  ph
  hold off;
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to L^2 norm}');
  xlabel('{\bf Mesh width}');
  ylabel('{\bf Discretization error}');
  
   add_Slope(gca,'SouthEast',1);
   add_Slope(gca,'SouthEast',2);
   add_Slope(gca,'SouthEast',3);
  clear all;