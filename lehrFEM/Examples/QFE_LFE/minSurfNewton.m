% Run script for minimum surface problem with the use of Newton's method

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

% Initialize constants

TOL = 1E-5;    % Error tolerence 
MAXIT = 100;    % Maximum number of iterations
NREFS = 2;      % Number of red refinement steps
PlotFreq = 3;  % Iterations between plots
GD_HANDLE = @(x,varargin)abs(acos(x(:,1))); % Boundary data
% Initial guess for square mesh
U0 = @(x,varargin)abs(acos(x(:,1))) + 10*sin(pi*x(:,1)/4).*sin(pi*x(:,2)/4);

% Initial guess for circular mesh
%U0 = @(x,varargin)x(:,1).*sin(x(:,2)) + sin(x(:,1)).^2+sin(x(:,2)).^2; % Initial guess of solution

DHANDLE = @dist_circ;
DAMPED = 1;     % Damped iteration / Armijo's rule
SMIN = 0.1;     % Minium value of s in damped iteration

% Initialize mesh


Mesh = load_Mesh('Coord_Circ.dat','Elem_Circ.dat'); % Remember to change input to refine_REG
%Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
%Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;
for i = 1:NREFS,
  Mesh = refine_REG(Mesh,DHANDLE,[0 0],1);
end
Mesh = add_Edge2Elem(Mesh);

% Interpolate initial data

M = assemMat_LFE(Mesh,@MASS_LFE);
L = assemLoad_LFE(Mesh,P3O3(),U0);
U_old = M\L;
plot_LFE(U_old,Mesh);
%colorbar;

err = 2*TOL;
it = 0;
while err > TOL,
    it = it + 1;
    if it > MAXIT,
        error('Too many iterations');
    end
    
    % Assemble matrices and load vector
    
    A = assemMAT_minSurf(U_old,Mesh,@STIMA_minSurf);
    B = assemMAT_minSurf(U_old,Mesh,@STIMA_minSurfDer);
    L = -A*U_old;
    
    % Solve the linear system
   
    [Delta,FreeDofs] = assemDir_LFE(Mesh,-1,@(x,varargin)0);
    Delta(FreeDofs) = B(FreeDofs,FreeDofs)\L(FreeDofs);
     
    err = norm(Delta(FreeDofs));
    
    % Damping strategy (Armijo's rule)
    
    if(DAMPED)
      
      s = 1;
      nrm = inf;
      while(s > SMIN && nrm > (1-s)*err)
        
        U_new = U_old + s*Delta;
        
        A = assemMat_minSurf(U_new,Mesh,@STIMA_minSurf);
        L = -A*U_new;

        [z,FreeDofs] = assemDir_LFE(Mesh,-1,@(x,varargin)0);
        z(FreeDofs) = B(FreeDofs,FreeDofs)\L(FreeDofs);
      
        nrm = norm(z);
        s = s/2;

      end
    
    else
        
      U_new = U_old+Delta;
        
    end
    
    U_old = U_new;
    
    % Plots every PlotFreq iterations
   
    if rem(it,PlotFreq) == 0,
       plot_LFE(U_new,Mesh);
       colorbar;
       title('{\bf Approximate solution}');
       xlabel(['{\bf ' int2str(it) ' iterations}']);
       drawnow();
              
    end
    
end
plot_LFEfancy(U_new,Mesh);