function varargout = main_psor(Shape,varargin)

% MAIN_PSOR Error analysis of the projected sor solver for LCP problems
%
%   MAIN_PSOR(SHAPE) Plots the error of the projected sor algorithm when
%   solving the obstackle problem with 5 red mesh refinements, 1e-10 as
%   iteration stopping criterion, uses MAXIT = 10000 as the maximum
%   number of iterations and sets OMEGA = 2/(1+sqrt(2)) as the 
%
%   SHAPE is the shape of the domain. Possible values are 'Sqr' for square
%   domain, 'Circ' for circular domain and 'LShap' for L shaped domain
%
%   MAIN_PSOR(SHAPE,NREFS) also let you specify the number of red
%   refinements of the mesh
%
%   MAIN_PSOR(SHAPE,NREFS,TOL) also specify the iteration stopping
%   criterion TOL
%
%   MAIN_PSOR(SHAPE,NREFS,TOL,MAXIT) specifies the maximum number of
%   iterations MAXIT
%
%   MAIN_PSOR(SHAPE,NREFS,TOL,MAXIT,OMEGA) let you specify the sor
%   extrapolation factor OMEGA
%
%   Example:
%
%   main_psor('Sqr');    

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

switch nargin
    case 1
        NREFS = 5;              % Number of red refinments
        TOL   = 1e-10;          % Stopping criterion
        MAXIT = 10000;          % Maximum number of iterations
        omega = 2/(1+sqrt(2));  % SOR extrapolation factor
    case 2
        NREFS = varargin{1};
        TOL   = 1e-10;      
        MAXIT = 10000;
        omega = 2/(1+sqrt(2));
    case 3
        NREFS = varargin{1};
        TOL   = varargin{2};
        MAXIT = 10000;
        omega = 2/(1+sqrt(2));
    case 4
        NREFS = varargin{1};
        TOL   = varargin{2};
        MAXIT = varargin{3};
        omega = 2/(1+sqrt(2));
    case 5
        NREFS = varargin{1};
        TOL   = varargin{2};
        MAXIT = varargin{3};
        omega = varargin{4};
end
    
% TOL = 1e-10;
% MAXIT = 5000;
% omega = 2/(1+sqrt(2));
% NREFS = 5;
GD_HANDLE = @(x,varargin)0;

switch Shape
    case 'Circ'
        % For circular mesh
        F_HANDLE = @(x,varargin)1;
        Obst = @(x,varargin)1/8+1/2*sqrt(x(:,1).^2+x(:,2).^2)-3/4*(x(:,1).^2+x(:,2).^2);
        Mesh = load_Mesh('Coord_Circ.dat','Elem_Circ.dat');
        type = 'circular';
    case 'Sqr'
        % For square mesh
        F_HANDLE = @(x,varargin)1;
        Obst = @(x,varargin)5*(1-x(:,1).^2).*(1-x(:,2).^2) + sin(pi*x(:,1)) .* sin(pi*x(:,2));
        Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
        type = 'square';
    case 'LShap'
        % For L-shaped mesh
        F_HANDLE = @(x,varargin)-5;
        Obst = @(x,varargin)5*(1-x(:,1).^2).*(1-x(:,2).^2).*x(:,1).*x(:,2);
        Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
        type = 'LShaped';
end



Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

for i = 1:NREFS,
    switch Shape
        case 'Circ'
            Mesh = refine_REG(Mesh,@dist_circ,[0 0],1);
        otherwise
            Mesh = refine_REG(Mesh);
    end
end
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

% Calculate solution

[x conv ERR it] = psor(A(FreeDofs,FreeDofs),b,x,omega,TOL,MAXIT);

if ~conv,
    error('The system does not converge within the specified number of iterations')
end

tol = TOL * ones(size(ERR,2),1);
fig = figure('Name', 'CG solver');
plot(ERR)
hold on
plot(tol,'r-');
hold off
set(gca,'YScale','log')
title(['{\bf Projected SOR solver for ', type, ' mesh}']);
xlabel('{\bf Iteration number}');
ylabel('{\bf Residual}');
name = ['psor_'  type  '.eps'];
print('-depsc', name);