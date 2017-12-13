% error_Obstacle2D_cg discretisation error analysis
%   Disretization error analysis of the obstacle problem using the
%   conjugate gradient method on a circular mesh with radius 1 and center
%   at the origin. The exact solution is given by
%   U = 1/4*(1-x(:,1).^2-x(:,2).^2) and the obstacle is 
%   Obst = 1/8+1/2*sqrt(x(:,1).^2+x(:,2).^2)-3/4*(x(:,1).^2+x(:,2).^2)
%   The plot of the error against mesh width is saved as error_cg.eps and
%   the plot of the error against the degrees of freedom is saved as
%   error_cg_DF.eps

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland


% Initialize constants

close all
clc

N = 5;
TOL = 1e-8;
MAXIT = 1e3;
GD_HANDLE = @(x,varargin)0;
F_HANDLE = @(x,varargin)1;
Obst = @(x,varargin)1/8+1/2*sqrt(x(:,1).^2+x(:,2).^2)-3/4*(x(:,1).^2+x(:,2).^2);
U_true = @(x,varargin)1/4*(1-x(:,1).^2-x(:,2).^2);
U_grad = @(x,varargin)-1/2*[x(:,1) x(:,2)];

% Prealocate memory

L2Err = zeros(N,1);
H1SErr = zeros(N,1);
h = zeros(N,1);
DF = zeros(N,1);

% Initialize mesh

Mesh = load_Mesh('Coord_Circ.dat','Elem_Circ.dat');
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;
for j = 1:N,
    Mesh = refine_REG(Mesh,@dist_circ,[0 0],1);
    Mesh = add_Edge2Elem(Mesh);
    
    % Assemble stiffness matrix and load vector

    A = assemMat_LFE(Mesh, @STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    
    % Assign obstacle
    
    nCoord = size(Mesh.Coordinates,1);
    psi = zeros(nCoord,1);
    for i = 1:nCoord,
        coord = Mesh.Coordinates(i,:);
        psi(i) = Obst(coord);
    end
    
    % Incorporate boundary data (Dirichlet)

    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    
    % Shift problem

    x = U(FreeDofs) - psi(FreeDofs);
    b = L(FreeDofs) - A(FreeDofs,FreeDofs)*psi(FreeDofs);
    
    % Calculate solution (with conjugate gradient method)

    [x conv] = cg(A(FreeDofs,FreeDofs),b,x,TOL,MAXIT);
    

    if ~conv,
        warning('The system does not converge within the specified number of iterations')
    end
    
    U(FreeDofs) = x + psi(FreeDofs);
    
    % Calculate discretisation error
    
    L2Err(j) = L2Err_LFE(Mesh,U,P7O6(),U_true);
    H1SErr(j) = H1SErr_LFE(Mesh,U,P7O6(),U_grad);
    h(j) = get_MeshWidth(Mesh);
    DF(j) = size(Mesh.Coordinates,1);
    
end

  % Plot out discretization error against h mesh width
  
  figure('Name','Discretization error against mesh width');
  plot(h,L2Err,'r-+', ...
       h,H1SErr,'b-+');
  legend('L2 error', ...
         'H1S error', ...
         'Location','SouthEast');
  title('{\bf Discretization errors}');
  xlabel('{\bf Mesh width [logarithmic]}');
  ylabel('{\bf Discretization error [logarithmic]}');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  legend('L2 norm','Location','SouthEast','H1 Semi-norm','Location','SouthEast');
  
  p = polyfit(log(h),log(H1SErr),1);
  add_Slope(gca,'NorthWest',p(1));
  p = polyfit(log(h),log(L2Err),1);
  add_Slope(gca,'SouthEast',p(1));
  
  print -depsc error_cg.eps
  
  % Plot out discretization error against DF degrees of freedom
  
  figure('Name','Discretization error against degrees of freedom');
  plot(DF,L2Err,'r-+', ...
       DF,H1SErr,'b-+');
  legend('L2 error', ...
         'H1S error', ...
         'Location','SouthEast');
  title('{\bf Discretization errors}');
  xlabel('{\bf Degrees of freedom [logarithmic]}');
  ylabel('{\bf Discretization error [logarithmic]}');
  set(gca,'XScale','log','YScale','log');
  legend('L2 norm','Location','SouthEast','H1 Semi-norm','Location','SouthEast');
  
  p = polyfit(log(DF),log(H1SErr),1);
  add_Slope(gca,'NorthEast',p(1));
  p = polyfit(log(DF),log(L2Err),1);
  add_Slope(gca,'SouthWest',p(1));
  
  print -depsc error_cg_DF.eps