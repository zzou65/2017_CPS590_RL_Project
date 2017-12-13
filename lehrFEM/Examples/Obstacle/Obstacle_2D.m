% Obstackle 2D
% Finite elements

% Initialize constants

TOL = 1e-5;
MAXIT = 100;
NREFS = 4;
isw = 1; % 1=PSOR, 2 = CG
Omega = 2/(1+sqrt(2));
GS_HANDLE = @(x,varargin)0;
F_HANDLE = @(x,varargin)-1;
Obst = @(x,varargin)5*cos(pi*x(:,1)).*(1-x(:,1).^2) .* (1-x(:,2).^2)-2;

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

% Assemble stiffness matrix, load vector and obstackle

A = assemMat_LFE(Mesh, @STIMA_Lapl_LFE);
L = assemLoad_LFE(Mesh,P706(),F_HANDLE);
psi = assemLoad_LFE(Mesh,P706(),Obst);

% Incorporate boundary data (Dirichlet)

[U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
L = L - A*U;

% Shift problem

x = U(FreeDofs) - psi(FreeDofs);
b = L(FreeDofs) - A(FreeDofs,FreeDofs)*psi(FreeDofs);

% Calculate solution

if isw==1,
    x = psor(A(FreeDofs,FreeDofs),b,x,Omega,TOL,MAXIT);
elseif isw==2,
    x = cg(A(FreeDofs,FreeDofs),b,x,TOL,MAXIT);
end

U(FreeDofs) = x + psi(FreeDofs);

plot_LFE(U,Mesh)