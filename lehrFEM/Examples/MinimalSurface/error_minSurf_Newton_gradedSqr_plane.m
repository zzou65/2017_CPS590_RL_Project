% error_minSurf_Newton_gradedSqr
%   Error analysis of the minimal surface problem on a graded square mesh
%   with with corner coordinates (-1,-1),(-1,1),(1,1),(-1,1) and grading
%   factor 2.5 towards the singularity points
%   Exact solution is given by U = x+y
%   The discretization error plot is saved as error_Newton_graded_plane.eps
%   The plot of the error against mesh width is saved as
%   error_Newton_graded_plane.eps the error against degrees of freedom is
%   saved as error_Newton_graded_plane_DF.eps

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland


clc
close all
% Initialize constants

TOL = 1E-12;    % Error tolerence 
MAXIT = inf;    % Maximum number of iterations
NREFS = 5;      % Number of different meshes

U_exact = @(x,varargin)x(:,1) + x(:,2);                       % Exact solution
GD_HANDLE = U_exact;                                                            % Boundary data
U0 = @(x,varargin)U_exact(x,varargin{:}) + (x(:,1).^2.*exp(x(:,2))).*sin(pi*x(:,1)).*sin(pi*x(:,2));   % Initial guess
U_grad = @(x,varargin)ones(size(x));

DAMPED = 1;    % Damped iteration / Armijo's rule
SMIN = 1e-3;   % Minium value of s in damped iteration

% Preallocate memory

L2Err = zeros(NREFS,1);
H1SErr = zeros(NREFS,1);
H = zeros(NREFS,1);
DF = zeros(NREFS,1);

for i = 1:NREFS,    
    % Initialize mesh
    
    Cname = ['Coord_Sqr_' int2str(i) '.dat'];
    Ename = ['Elem_Sqr_' int2str(i) '.dat'];
    Mesh = load_Mesh(Cname,Ename);
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh = add_Edge2Elem(Mesh);
 
    % Interpolate initial data
    
    M = assemMat_LFE(Mesh,@MASS_LFE);
    L = assemLoad_LFE(Mesh,P3O3(),U0);
    
    U_old = M\L;
    
    err = 2*TOL;
    it = 0;
    
    while err>TOL,
        it = it + 1;
        if it>MAXIT,
            error('Tpp many iterations');
        end
        
        % Assemble matrices and load vector
        
        A = assemMat_minSurf(U_old,Mesh,@STIMA_minSurf);
        B = assemMat_minSurf(U_old,Mesh,@STIMA_minSurfDer);
        L = -A*U_old;
        
        % Solve the linear system
        
        [Delta,FreeDofs] = assemDir_LFE(Mesh,-1,@(x,varargin)0);
        Delta(FreeDofs) = B(FreeDofs,FreeDofs)\L(FreeDofs);
        
        err = norm(Delta(FreeDofs));
        
        % Damping stategy (Armijo's rule)
        
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
            U_new = U_old + Delta;
        end
        
        U_old = U_new;
    end
    for j = 1:size(Mesh.Coordinates,1),
        U_ex(j) = U_exact(Mesh.Coordinates(j,:));
    end
    
    L2Err(i) = L2Err_LFE(Mesh,U_new,P7O6(),U_exact);
    H1SErr(i) = H1SErr_LFE(Mesh,U_new,P7O6(),U_grad);
    H(i) = get_MeshWidth(Mesh);
    DF(i) = size(Mesh.Coordinates,1);
end


% Plot of error against mesh width

figure('Name','Discretization error');
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

print -depsc error_Newton_graded_plane.eps

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

print -depsc error_Newton_graded_plane_DF.eps