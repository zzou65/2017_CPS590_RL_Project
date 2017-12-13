% error_minSurf_Newton

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

% Initialize constants
clc
close all
K = 50;
TOL = 1E-12;    % Error tolerence 
MAXIT = inf;    % Maximum number of iterations
NREFS = 4;      % Number of red refinement steps
GD_HANDLE = @(x,varargin)((cosh(x(:,2))).^2 - x(:,1).^2).^(1/2); % Boundary data
% Initial guess for square mesh
func = @(x,varargin)5*sin(x(:,1)*pi)+10*sin(x(:,2)*pi/5);
U_1 = @(x,varargin)K/sinh(pi/2) * sin(pi/2*x(:,1)) .* sinh(pi/2*(1-x(:,2)));
%U0 = @(x,varargin)U_1([x(1) 0]) + x(2) .* (U_1([x(1) 1]) - U_1([x(1) 0]));
%U0 = @(x,varargin)sin(pi*x(:,1)) .* (1-x(:,2));
U0 = @(x,varargin)(1-x(:,1).^2) + x(:,2) .* (sqrt(cosh(1).^2 - x(:,1).^2) - sqrt(1 - x(:,1).^2));
%U0 = @(x,varargin)func(x) .* sin(pi*x(:,1)) .* sin(pi*x(:,2)) + ((cosh(x(:,2))).^2 - x(:,1).^2).^(1/2); % Initial guess
U_exact = @(x,varargin)((cosh(x(:,2))).^2 - x(:,1).^2).^(1/2);

DAMPED = 1;     % Damped iteration / Armijo's rule
SMIN = 0.1;     % Minium value of s in damped iteration

% Initialize mesh


Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;
for i = 1:NREFS,
  Mesh = refine_REG(Mesh);
end
Mesh = add_Edge2Elem(Mesh);
% Interpolate initial data

M = assemMat_LFE(Mesh,@MASS_LFE);
L = assemLoad_LFE(Mesh,P3O3(),U0);

U_old = M\L;
%plot_LFE(U_old,Mesh);
%colorbar;

err = 2*TOL;
it = 0;
while err > TOL,
    it = it + 1;
    if it > MAXIT,
        error('Too many iterations');
    end
    
    % Assemble matrices and load vector
    
    A = assemMat_minSurf(U_old,Mesh,@STIMA_minSurf);
    B = assemMat_minSurf(U_old,Mesh,@STIMA_minSurfDer);
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
    
end

U_ex = zeros(size(Mesh.Coordinates,1),1);
nCoord = size(Mesh.Coordinates,1);
for i = 1:nCoord,
    coord = Mesh.Coordinates(i,:);
    U_ex(i) = U_exact(coord);
end



L2Err = L2Err_LFE(Mesh,U_new,P7O6(),U_exact)