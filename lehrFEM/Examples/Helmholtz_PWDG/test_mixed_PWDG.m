function [] = test_mixed_PWDG()
%TEST_MIXED_PWDG solve sample problem with mixed PWDG
%   Detailed explanation goes here

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Define problem parameters
omega = 1.5*pi;
d = [1 0.5];
d = d/norm(d);
lambda = 0.5*omega;
psi = 0;
u_ex = @(x,varargin) sin(lambda*x*d'+psi);
du_ex = @(x,varargin) lambda*cos(lambda*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);
gI = @(x,n,varargin) lambda*(d*n')*cos(lambda*x*d'+psi) + i*omega*sin(lambda*x*d'+psi);
f = @(x,varargin) (lambda^2-omega^2)*sin(lambda*x*d'+psi);

% Define discretization parameters
nref = 4;
ndir = 5;

qr1 = gauleg(0,1,12);
qr2 = Duffy(TProd(qr1));
% qr2 = TProd(qr1);

% Construct mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
Mesh.Elements = [1 2 3; 1 4 3];
% Mesh.Elements = [1 2 3 4];
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
% Mesh = set_Data_PWDG(Mesh,'Dir',[d;-d],'Omega',omega,...
%   'a',@(w,h,varargin) 2/(w*h),'b',@(w,h,varargin) 0.01*w*h,'d',0.75);
Mesh = set_Data_PWDG(Mesh,'nDofs',ndir,'MakeDir','uniform','Omega',omega,...
  'a',@(w,h,varargin) 2/(w*h),'b',@(w,h,varargin) 0.01*w*h);

% Assemble L2 inner product
[I_mass_inn,J_mass_inn,M_inn] = assemMat_Inn_PDG2(Mesh,@MASS_Inn_PWDG,omega);
[I_mass_bnd,J_mass_bnd,M_bnd] = assemMat_Bnd_PDG2(Mesh,[],@MASS_Bnd_PWDG,omega);
M = sparse([I_mass_inn;I_mass_bnd],[J_mass_inn;J_mass_bnd],[M_inn;M_bnd]);

% Assemble stiffness matrix
[I_Vol,J_Vol,A_Vol] = assemMat_Vol_PDG2_vec(Mesh,3,@STIMA_Helm_Vol_PWDG_dual,omega,M);
[I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2_vec(Mesh,3,@STIMA_Helm_Inn_PWDG_dual);
[I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2_vec(Mesh,[],3,@STIMA_Helm_Imp_Bnd_PWDG_dual);
A_m = sparse([I_Vol;I_Inn;I_Bnd],[J_Vol;J_Inn;J_Bnd],[A_Vol;A_Inn;A_Bnd]);

% Assemble load vector
b_g = assemLoad_Bnd_PDG2_vec(Mesh,[],3,@LOAD_Imp_Bnd_PWDG_dual,qr1,omega,gI);
b_f = assemLoad_Vol_PDG2_vec(Mesh,3,@LOAD_Vol_PWDG_dual,qr2,omega,f);
b_m = b_f + b_g;

% Solve equation
U = A_m\b_m;

% Extract solution
u_m = U(1:end/3);

% Compute and print errors
err0_m = L2Err_PWDG(Mesh,u_m,qr2,omega,u_ex);
err1_m = EnergyErr_PWDG(Mesh,u_m,qr2,omega,u_ex,du_ex,0);
fprintf('Mixed PWDG Error: L2: %g    En: %g\n',err0_m,err1_m);


% Set sigma to grad(u)/(i*omega) and solve restricted equation
%   A_r is i*omega times stiffness matrix of primal method
Dir = vertcat(Mesh.ElemData.Dir);
n = size(Dir,1);
P = [speye(n);spdiags(Dir(:,1),0,n,n);spdiags(Dir(:,2),0,n,n)];
A_r = P'*A_m*P;
b_r = P'*b_m;
u_r = A_r\b_r;

% Compute and print errors
err0_r = L2Err_PWDG(Mesh,u_r,qr2,omega,u_ex);
err1_r = EnergyErr_PWDG(Mesh,u_r,qr2,omega,u_ex,du_ex,0);
fprintf('Primal PWDG Error: L2: %g    En: %g\n',err0_r,err1_r);

% Plot solution
plot_PWDG(u_m,Mesh,omega,max(0,5-nref));

% Add new line
fprintf('\n');

return