% Computes the discretization error measured in the H1 semi-norm and the
% L2norm at the terminal time using the theta scheme to approximate the
% solution to the parabolic problem.
%
% This script file generates the .eps figures:
%  Parabolic_L2MaxErr.eps,
%  Parabolic_L2MeanErr.eps
%  Parabolic_H1SMeanErr.eps

% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NMESHREFS = 7;                                                       % Number of red mesh refinements
  NSTEPREFS = 11;                                                      % Number of time step refinements
  T = 21.3;                                                            % Final time
  THETA = 1;                                                           % Coefficient of theta scheme
  GD = @(x,varargin)zeros(size(x,1),1);                                % Dirichlet boundary data
  F = @(x,Elemflag,t,varargin)2*pi*(4*pi*cos(2*pi*t)-sin(2*pi*t)) ...
                              *sin(2*pi*x(:,1)).*sin(2*pi*x(:,2));     % Right hand side source term        
  U0 = @(x,varargin)(2*sin(2*pi*x(:,1)).*sin(2*pi*x(:,2)));            % Initial data
  UEX = @(x,t,varargin)sin(2*pi*x(:,1)).*sin(2*pi*x(:,2)) ...
                       .*(cos(2*pi*t)+exp(-8*pi^2*t));                 % Exact solution
  GRAD_UEX = @(x,t,varargin)2*pi*(cos(2*pi*t)+exp(-8*pi^2*t))*...      % Gradient of exact solution 
                            [cos(2*pi*x(:,1)).*sin(2*pi*x(:,2)) ...
                             sin(2*pi*x(:,1)).*cos(2*pi*x(:,2))];
  
  % Compute convergence rate
  
  fprintf('### Convergence rates for parabolic problems ');
  if(THETA == 0)
    scheme = 'Expl. Euler';  
  elseif(THETA == 1/2)
    scheme = 'Crank-Nicolson';  
  elseif(THETA == 1)
    scheme = 'Imp. Euler';  
  else
    scheme = ['theta = ' num2str(theta,'%.1f')]; 
  end
  fprintf(['(' scheme ')\n\n']); 
  
  % Initialize variables
  
  h = zeros(1,NMESHREFS);
  dt = zeros(1,NSTEPREFS);
  L2MaxErr = zeros(NSTEPREFS,NMESHREFS);
  L2MeanErr = zeros(NSTEPREFS,NMESHREFS); 
  H1SMeanErr = zeros(NSTEPREFS,NMESHREFS);
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat'); 
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  
  for j1 = 1:NMESHREFS
        
    % Do one red refinement step
    
    Mesh = refine_REG(Mesh);
      
    % Precompute matrices
  
    [I1,J1,M] = assemMat_LFE(Mesh,@MASS_LFE);
    [I2,J2,A] = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    
    % Compute initial data
  
    u_old = assemLoad_LFE(Mesh,P1O2(),U0);
    u_old = sparse(I1,J1,M)\u_old;
    L_old = assemLoad_LFE(Mesh,P1O2(),F,0);
    
    nSteps = 4;
    for j2 = 1:NSTEPREFS      
      fprintf('   Number of red mesh refinements                   :  %8d\n',j1);
      fprintf('   Number of time step refinements                  :  %8d\n',nSteps);
      t = cputime();
      
      % Compute system matrices
      
      dt(j2) = T/nSteps;
      S1 = sparse([I1; I2],[J1; J2],[M; dt(j2)*THETA*A]);
      S2 = sparse([I1; I2],[J1; J2],[M; -dt(j2)*(1-THETA)*A]);
   
      % Start theta scheme
 
      L_old = assemLoad_LFE(Mesh,P1O2(),F,0);
      fprintf('   Running theta scheme                             :  ');
      per = 0;
      progress_bar(per);
      H1SErr = zeros(1,nSteps);
      L2Err = zeros(1,nSteps);
      for i = 1:nSteps
        if(per < floor(100*i/nSteps))
          per = floor(100*i/nSteps);
          progress_bar(per);
        end
      
	    ti = i*dt(j2);
    
        % Assemble load vector 
    
        L_new = assemLoad_LFE(Mesh,P1O2(),F,i*dt(j2));
    
        % Incorporate Dirichlet boundary data  
      
        [u_new,FreeDofs] = assemDir_LFE(Mesh,-1,GD,i*dt(j2));
    
        % Solve the linear system
    
        rhs = S2*u_old + dt(j2)*THETA*L_new+dt(j2)*(1-THETA)*L_old;
        rhs = rhs - S1*u_new;
        u_new(FreeDofs) = S1(FreeDofs,FreeDofs)\rhs(FreeDofs);
    
	    % Compute discretization errors 
    
	    H1SErr(i) = H1SErr_LFE(Mesh,u_new,P1O2(),GRAD_UEX,i*dt(j2));
	    L2Err(i)  = L2Err_LFE(Mesh,u_new,P1O2(),UEX,i*dt(j2));
        	
        u_old = u_new;
        L_old = L_new;  
      end
      fprintf('    done\n');
           
      % Update time steps (half size of timestep)

      nSteps = 2*nSteps;
      
      % Collect overall errors
      
      H1SMeanErr(j2,j1) = sqrt(dt(j2)*sum(H1SErr.^2));
      L2MaxErr(j2,j1) = max(L2Err);
      L2MeanErr(j2,j1)  = sqrt(dt(j2)*sum(L2Err.^2));
            
      fprintf('   Overall runtime at current refinement level [m]  :  %8.2f\n\n',(cputime()-t)/60);
    end
    
    % Compute mesh with and degrees of freedom
    
    h(j1) = get_MeshWidth(Mesh);
    
  end
  
  % Generate data files and figures
  
  [h,dt] = meshgrid(h,dt);
  
  fig = figure('Name','Discretization error');
  surf(h,dt,L2MaxErr, ...
       'EdgeColor','r', ...
       'FaceColor','w');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Time step [log]}');
  title('{\bf Maximal discretization error wrt. L^2 norm}');
  set(gca,'XScale','log','YScale','log','ZScale','log');
  print('-depsc','Parabolic_L2MaxErr.eps');
 
  fig = figure('Name','Discretization error');
  surf(h,dt,L2MeanErr, ...
       'EdgeColor','r', ...
       'FaceColor','w');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Time step [log]}');
  title('{\bf Mean discretization error wrt. L^2 norm}');
  set(gca,'XScale','log','YScale','log','ZScale','log');
  print('-depsc','Parabolic_L2MeanErr.eps');
  
  fig = figure('Name','Discretization error');
  surf(h,dt,H1SMeanErr, ...
       'EdgeColor','r', ...
       'FaceColor','w');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Time step [log]}');
  title('{\bf Mean discretization error wrt. H^1 semi-norm}');
  set(gca,'XScale','log','YScale','log','ZScale','log');
  print('-depsc','Parabolic_H1SMeanErr.eps');
  
  % Save data to file
  
  save('Parabolic.mat','h','dt','L2MaxErr','L2MeanErr','H1SMeanErr');
  
  % Clear memory
  
  clear all;
  