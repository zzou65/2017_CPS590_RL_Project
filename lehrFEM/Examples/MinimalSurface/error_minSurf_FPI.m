% error_minSurf_FPI
%   Error analysis of the minimal surface problem using fixed point
%   iteration on a uniform square mesh with corner coordinates 
%   (0,0),(1,0),(1,1),(0,1). The plot of the error against mesh width is
%   saved as error_FPI.eps and the plot of the error against degrees of
%   freedom is saved as error_FPI_DF.eps
%   Exact solution is given by U = sqrt(cosh(y).^2 - x.^2) and the gradient
%   of U is discontinious

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

clc
close all
% Initialize constants

TOL = 1e-12;
MAXIT = inf;
NREFS = 6;

U_exact = @(x,varargin)sqrt(cosh(x(:,2)).^2 - x(:,1).^2);                       % Exact solution
GD_HANDLE = U_exact;                                                            % Boundary data
U0 = @(x,varargin)U_exact(x,varargin{:}) + .5.*sin(pi*x(:,1)).*sin(pi*x(:,2));  % Initial guess
U_grad = @Ugrad;                                                                % Gradient of soulution


% Initialize mesh

Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

% Preallocate memory

L2Err = zeros(NREFS,1);
H1SErr = zeros(NREFS,1);
H = zeros(NREFS,1);
DF = zeros(NREFS,1);

for i = 1:NREFS,
    Mesh = refine_REG(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    
    % Interpolate initial data

    M = assemMat_LFE(Mesh,@MASS_LFE);
    L = assemLoad_LFE(Mesh,P3O3(),U0);
    U = M\L;
    U_0 = U;
    
    err = 2 * TOL;
    it = 0;
    
    while err > TOL,
       it = it + 1;
       if it > MAXIT,
           error('Too many iterations');
       end
       U_old = U; % Keep this for termination criteria

       A = assemMat_minSurf(U,Mesh,@STIMA_minSurf);
       L = zeros(size(Mesh.Coordinates,1),1);
       
       % Incorporate Dirichlet boundary data
       
       [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
       L = L - A*U;
       
       % Solve the linear system

       U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
       
       % Check if we are done

       err = norm(U(FreeDofs)-U_old(FreeDofs))/norm(U_0(FreeDofs));
    end
    L2Err(i) = L2Err_LFE(Mesh,U,P7O6(),U_exact);
    H1SErr(i) = H1SErr_LFE(Mesh,U,P7O6(),U_grad);
    H(i) = get_MeshWidth(Mesh);
    DF(i) = size(Mesh.Coordinates,1);
end    

% Plot of error against mesh width

figure('Name','Discretization error against mesh width');
plot(H,L2Err,'r-+', ...
     H,H1SErr,'b-+')
legend('L2 error', ...
       'H1S error', ...
       'Location', 'SouthEast');
title('{\bf Discretization errors}');
xlabel('{\bf Mesh width [logarithmic]}');
ylabel('{\bf Discretization error [logarithmic]}');
set(gca, ...
    'XDir','reverse', ...
    'XScale','log', ...
    'YScale','log');
legend('L2 norm','Location', 'SouthEast','H1 Semi-norm','Location','SouthEast');

p = polyfit(log(H),log(H1SErr),1);
add_Slope(gca,'NorthWest',p(1));
p = polyfit(log(H),log(L2Err),1);
add_Slope(gca,'SouthEast',p(1));

print -depsc error_FPI.eps

% Plot of error against degrees of freedom


figure('Name','Discretization error against degrees of freedom');
plot(DF,L2Err,'r-+', ...
     DF,H1SErr,'b-+')
legend('L2 error', ...
       'H1S error', ...
       'Location', 'SouthEast');
title('{\bf Discretization errors}');
xlabel('{\bf Degrees of freedom [logarithmic]}');
ylabel('{\bf Discretization error [logarithmic]}');
set(gca, ...
    'XScale','log', ...
    'YScale','log');
legend('L2 norm','Location', 'SouthEast','H1 Semi-norm','Location','SouthEast');

p = polyfit(log(DF),log(H1SErr),1);
add_Slope(gca,'NorthEast',p(1));
p = polyfit(log(DF),log(L2Err),1);
add_Slope(gca,'SouthWest',p(1));

print -depsc error_FPI_DF.eps