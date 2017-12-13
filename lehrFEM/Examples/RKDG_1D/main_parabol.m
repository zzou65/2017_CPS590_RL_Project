% Run script for linear advection problem.

% Copyright 2006-2007 Patrick Meury, Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland
  
  % Initialize constants
  X0 =-1;                 % Left end point of interval
  X1 = 1;                 % Right end point of interval
  NPOINTS = 11;           % Initial set of grid points
  NREFS = 5;              % Number of red refinement steps
  P=0;                    % ploynomial degree
  NSTEPS =300;             % Number of time steps
  AVIFILE = 'Shock.avi';  % Filename of .avi file           
  NFRAMES = 00;          % Number of movie frames
  T =40;                   % Final time
  s=0;                   % switch between: 1 SIP, -1 NIP, 0 IIP  
  alpha=2;               % stabilizaton term   
  
  U0 = @(x,t)x.^2-1;                            % Initial data
  UD = @(x,t)(x.^2-1)*cos(pi*t);                % Dirichlet boundary data
  G = @(x,t)-(x.^2-1)*pi*sin(pi*t)-2*cos(pi*t);  % Right hand side load data
  FLUX = @(x)zeros(size(x,1));                  % Flux function
  VHandle =@(x)zeros(size(x,1));
  NUMFLUX = @Upwind_LA;   % Numerical flux function
  ULIM = @limit_ad;       % Limiter for ansatz functions
  
  % Initialize TVD scheme (2-stage Runge-Kutta scheme)
%     
  NSTAGES = 2;
  A_RK = [1 0; 1/2 1/2];
  B_RK = [1 0; 0   1/2];
%   NSTAGES = 1;
%   A_RK = [1];
%   B_RK = [1];
%NSTAGES=3;
%    A_RK=RKDG_alpha(NSTAGES);
%    B_RK=RKDG_beta(NSTAGES);
  
%   % Open FIFO buffer
%   
%   buf = open(buffer());
%   
  % Initialize mesh
  
  Coordinates = X0 + (X1-X0)/(NPOINTS-1)*(0:(NPOINTS-1));

  % Initialize polynomial degrees
  
  p=P*ones(1,NPOINTS-1);
  dofs=sum(p+1);
  
  % quadrature rule and shape functions at quadrature point and boundary
  
  QuadRule=gauleg(-1,1,100);
  
  Shap=shap_Leg_1D(QuadRule.x,P);
  Grad_Shap=grad_shap_Leg_1D(QuadRule.x,P);
  Inn_shap=shap_Leg_1D([-1;1],P);
  Inn_gradshap=grad_shap_Leg_1D([-1;1],P);
  
%   Grad_Shap = grad_shap_DGLFE_1D(QuadRule.x);
%   Shap = shap_DGLFE_1D(QuadRule.x);
%   Inn_shap=shap_DGLFE_1D([-1;1]);
%   Inn_gradshap=grad_shap_DGLFE_1D([-1;1]);
  
  %Assemble mass matrix and extract diagonal
  M=assemMat_Vol_hpDG_1D(Coordinates,p,@MASS_Vol_hpDG_1D,QuadRule,Shap);
  M=spdiags(M,0);
  
  A_vol=assemMat_Vol_hpDG_1D(Coordinates,p,@STIMA_Vol_Lapl_hpDG_1D,QuadRule,Grad_Shap);
  A_Inn=assemMat_Inn_hpDG_1D(Coordinates,p,@STIMA_Inn_Lapl_hpDG_1D,Inn_shap,Inn_gradshap,s);
  A_InnPen=assemMat_Inn_hpDG_1D(Coordinates,p,@STIMA_InnPen_hpDG_1D,Inn_shap,alpha);
  A_Bnd=assemMat_Bnd_hpDG_1D(Coordinates,p,Inn_shap,Inn_gradshap,s,alpha);
   
  %S=(A_vol-A_Inn-A_InnPen);
  S=(A_vol-A_Inn-A_Bnd+A_InnPen);
  % Compute initial data
  
  L = assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,U0,0);
  
  U = L./M;
  CMin = min(U);
  CMax = max(U);
%  buf = push(buf,U);
  
  % Integrate ODE system
    
  per = 0;
  progress_bar(per)    
  dt = T/NSTEPS;
  tau = zeros(1,NSTAGES+1);
  L = zeros(size(U,1),NSTAGES+1);
  V = zeros(size(U,1),NSTAGES+1);
  PV = zeros(size(U,1),NSTAGES+1);
  for i = 1:NSTEPS 
      
    if(per < floor(100*i/NSTEPS))
      per = floor(100*i/NSTEPS);
      progress_bar(per);
    end 
      
    tau(1) = (i-1)*dt; 
    V(:,1) = U;
   
    Lvol=assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,G,tau(1));
    
   % PV(:,1) = ULIM(Mesh,V(:,1),ALPHA);
     D=assemDir_hpDG_1D(Coordinates,p,UD,Inn_shap,Inn_gradshap,s,alpha,tau(1));
    
%      Bvol=assemRKDG_update_Vol(Coordinates,U,p,shap,gradshap,flowhandle,QuadRule);
%      Binn=assemRKDG_update_Inn(Coordinates,U,p,lshap,rshap,numflowhandle,flowhandle);
%      Bbnd=assemRKDG_update_Inn(Coordinates,U,p,lshap,rshap,numflowhandle,flowhandle);
%     
    L(:,1) = Lvol-S*V(:,1);
    %L(:,1)=Lvol;
    for j = 2:(NSTAGES+1)
      tau(j) = 0;
      V(:,j) = 0;

      for k = 1:(j-1)
        tau(j) = tau(j) + A_RK(j-1,k)*tau(k) + dt*B_RK(j-1,k);
        V(:,j) = V(:,j) + A_RK(j-1,k)*V(:,k) + dt*B_RK(j-1,k)*(L(:,k)./M);
      end
      
      if(j < NSTAGES+1)
        Lvol=assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,G,tau(j));
        D=assemDir_hpDG_1D(Coordinates,p,UD,Inn_shap,Inn_gradshap,s,alpha,tau(j)); 

  %     Bvol=assemRKDG_update_Vol(Coordinates,u,p,shap,gradshap,flowhandle,QuadRule);
  %     Binn=assemRKDG_update_Inn(Coordinates,u,p,lshap,rshap,numflowhandle,flowhandle);
  %     Bbnd=assemRKDG_update_Inn(Coordinates,u,p,lshap,rshap,numflowhandle,flowhandle);
    
        L(:,j) = Lvol-S*V(:,j);  
        %L(:,j) = Lvol;
      end 
      
    end
         
    % Update solution
    
    U = V(:,NSTAGES+1);
    
    if (mod(i,50)==0)
     plot_hpDG_1D(Coordinates,p,U,@shap_Leg_1D);
     %plot_hpDG_1D(Coordinates,p,U,@shap_DGLFE_1D)
     set(gca, ...
        'XLim',[-1,1], ...
        'YLim',[-2,2]);
    end
%     if(rem(i,max(floor(NSTEPS/NFRAMES),1)) == 0)
%       buf = push(buf,U);
%       CMax = max(CMax,max(U));
%       CMin = min(CMin,min(U));  
%     end

  end
  hold off;
  clear('M','V','PV','L','Lvol','Bvol','Binn','Bbnd');
  
  % Generate auxiliary mesh
  
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  Coordinates = zeros(nCoordinates,2);
  Elements = zeros(nElements,3);
  V = zeros(nCoordinates,1);
  
  for i = 1:nElements
    Elements(i,:) = 3*(i-1)+[1 2 3];
    Coordinates(Elements(i,1),:) = Mesh.Coordinates(Mesh.Elements(i,1),:);
    Coordinates(Elements(i,2),:) = Mesh.Coordinates(Mesh.Elements(i,2),:);
    Coordinates(Elements(i,3),:) = Mesh.Coordinates(Mesh.Elements(i,3),:);
  end
  
  % Generate movie
  
  if(~isempty(AVIFILE))
    Mov = avifile(AVIFILE);
  end
  
  fig = figure('Name','DGFEm for parabolic problems');
  XMax = max(Mesh.Coordinates(:,1));
  XMin = min(Mesh.Coordinates(:,1));
  XLim = [XMin XMax] + 0.05*(XMax-XMin)*[-1 1];
  YMax = max(Mesh.Coordinates(:,2));
  YMin = min(Mesh.Coordinates(:,2));
  YLim = [YMin YMax] + 0.05*(YMax-YMin)*[-1 1];
  CLim = [CMin CMax] + 0.05*(CMax-CMin)*[-1 1];
  
  pos = 1;
  while(~isempty(buf))
     
    % Extract solution from FIFO buffer  
       
    [U,buf] = pop(buf);
     
    V(Elements(:,1)) = U(Elements(:,2))+U(Elements(:,3))-U(Elements(:,1));
    V(Elements(:,2)) = U(Elements(:,1))+U(Elements(:,3))-U(Elements(:,2));
    V(Elements(:,3)) = U(Elements(:,1))+U(Elements(:,2))-U(Elements(:,3));
    
    % Generate figure
     
    clf;
    patch('Faces', Elements, ...
          'Vertices', [Coordinates(:,1) Coordinates(:,2) V], ...
          'CData', V, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    set(gca,'XLim',XLim, ...
            'YLim',YLim, ...
            'CLim',CLim, ...
            'DataAspectRatio',[1 1 1]);
    colorbar;
    drawnow();
    
    % Add frame to movie
    
    if(~isempty(AVIFILE))
      try
        F = getframe();
        Mov = addframe(Mov,F);
      catch          
        fprintf('Invalid frame at position %d\n',pos);  
      end
    end
    
    pos = pos+1;
    
   end
   buf = close(buf);
    
   % Close the movie
  
  if(~isempty(AVIFILE))
    Mov = close(Mov);
  end
    
  % Clear memory
  
  clear all;