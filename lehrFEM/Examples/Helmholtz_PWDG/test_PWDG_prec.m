function [] = test_PWDG_prec()
%TEST_PWDG_PREC solve sample problem with preconditioned PWDG

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Define problem parameters
omega = pi;
d = [1 0.4];
d = d/norm(d);
lambda = 0.9*omega;
psi = 0;
u_ex = @(x,varargin) sin(lambda*x*d'+psi);
du_ex = @(x,varargin) lambda*cos(lambda*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);
gI = @(x,n,varargin) lambda*(d*n')*cos(lambda*x*d'+psi) + i*omega*sin(lambda*x*d'+psi);
% gD = @(x,varargin) sin(lambda*x*d'+psi);
f = @(x,varargin) (lambda^2-omega^2)*sin(lambda*x*d'+psi);

% Define discretization parameters
nref = 2;
ndir = 21;
phi = 2*pi*(0:1/ndir:1-1/ndir)';
dir = [cos(phi),sin(phi)];

qr1 = gauleg(0,1,20);
qr2 = TProd(gauleg(0,1,20));
qr2 = Duffy(qr2);          % for triangles

% Construct mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
Mesh.Elements = [1 2 3; 1 4 3];               % triangular mesh
% Mesh.Elements = [1 2 3 4];                    % quadrilateral mesh
Mesh = add_Edges(Mesh);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;
for j = 1:nref
  Mesh = refine_REG(Mesh);
end
Mesh = orient_Elems(Mesh);
Mesh = add_Edge2Elem(Mesh);
Mesh = add_DGData(Mesh);
Mesh = set_Data_PWDG(Mesh,'Dir',dir,'Omega',omega,...
  'a',@(w,h,varargin) 1/(w*h),'b',0,'d',1);
% plot_Mesh(Mesh,'petas');

%%% Assemble stiffness matrix and load vector

% Assemble stiffness matrix
[I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
[I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,-1,@STIMA_Helm_Imp_Bnd_PWDG,omega);
% [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,-1,@STIMA_Helm_Dir_Bnd_PWDG,omega);
A = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

% Assemble load vector
b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
% b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Dir_Bnd_PWDG,qr1,omega,gD);
b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,f);
b = b_f + b_g;

% Solve equation without preconditioning
u0 = A\b;

% Construct preconditioner
[PL,PR] = assemPrec_SVD_PDG(Mesh,A);

% Solve equation with preconditioning
w = (PL*A*PR)\(PL*b);
u = PR*w;

% Plot solution
plot_PWDG(u,Mesh,omega,max(0,5-nref));

% Assemble L2 inner product
[I_mass_inn,J_mass_inn,M_inn] = assemMat_Inn_PDG2(Mesh,@MASS_Inn_PWDG,omega);
[I_mass_bnd,J_mass_bnd,M_bnd] = assemMat_Bnd_PDG2(Mesh,[],@MASS_Bnd_PWDG,omega);
M = sparse([I_mass_inn;I_mass_bnd],[J_mass_inn;J_mass_bnd],[M_inn;M_bnd]);

% Calculate L2 projection of exact solution onto discretization space
b_u = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,u_ex);
u_p = M\b_u;

% Calculate L2 discretization error
err = L2Err_PWDG(Mesh,u,qr2,omega,u_ex);
err0 = L2Err_PWDG(Mesh,u0,qr2,omega,u_ex);
err_p = L2Err_PWDG(Mesh,u_p,qr2,omega,u_ex);
fprintf('L2 Error : %g (prec.)    %g (no prec.)    L2 Proj.Error : %g\n',err,err0,err_p);

% Assemble energy-norm inner product
[I_h1_inn,J_h1_inn,H_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega);
[I_h1_bnd,J_h1_bnd,H_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega);
B = sparse([I_h1_inn;I_h1_bnd],[J_h1_inn;J_h1_bnd],[H_inn;H_bnd]);

% Compute energy-norm projection of exact solution onto PWDG-space
b_u = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega,u_ex,du_ex);
u_p = B\b_u;

% Calculate energy norm discretization error
err = EnergyErr_PWDG(Mesh,u,qr2,omega,u_ex,du_ex,0);
err0 = EnergyErr_PWDG(Mesh,u0,qr2,omega,u_ex,du_ex,0);
err_p = EnergyErr_PWDG(Mesh,u_p,qr2,omega,u_ex,du_ex,0);
fprintf('En.Error : %g (prec.)    %g (no prec.)    En.Proj.Error : %g\n\n',err,err0,err_p);

return