function [] = test_local()
%TEST_LOCAL Summary of this function goes here
%   Detailed explanation goes here

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Define problem parameters
omega = 1e-2;
d1 = [-1 0.001];
d1 = d1/norm(d1);
c1 = 1;
d2 = [1 1];
d2 = d2/norm(d2);
c2 = 0.35;
psi = 0;
u_ex = @(x,varargin) c1*exp(i*(omega*x*d1'+psi))+c2*exp(i*(omega*x*d2'+psi));
du_ex = @(x,varargin) c1*i*omega*exp(i*(omega*x*d1([1 1],:)'+psi)).*d1(ones(size(x,1),1),:)...
  + c2*i*omega*exp(i*(omega*x*d2([1 1],:)'+psi)).*d2(ones(size(x,1),1),:);
g = @(x,n,varargin) c1*i*omega*(1+d1*n')*exp(i*(omega*x*d1'+psi)) ...
  + c2*i*omega*(1+d2*n')*exp(i*(omega*x*d2'+psi));

% Define discretization parameters
ndir = 11;
phi = 2*pi*(0:1/ndir:1-1/ndir)';
dir = [cos(phi),sin(phi)];
ndir = ndir+2;
dir = [d1;d2;dir];
% ndir = 3;
% phi = pi*([0.2;1-1/42;1+1/19]);
% dir = [cos(phi),sin(phi)];

qr1 = gauleg(0,1,25);
qr2 = TProd(qr1);

% Construct mesh
Mesh.Coordinates = [0 0; 1 0; 1 1; 0 1];
Mesh.Elements = [1 2 3 4];
Mesh = add_Edges(Mesh);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;
Mesh = orient_Elems(Mesh);
Mesh = add_Edge2Elem(Mesh);
Mesh = add_DGData(Mesh);
Mesh = set_Data_PWDG(Mesh,'Dir',dir,'Omega',omega);
disp(Mesh.ElemData(1).Dir);

% Assemble stiffness matrix
[I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
[I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,-1,@STIMA_Helm_Imp_Bnd_PWDG,omega);
A = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

% Assemble load vector
b = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,g);

% Solve equation
u = A\b;

% % Assemble L2 inner product
% [I_mass_diag,J_mass_diag,M_diag] = assemMat_Vol_PDG2(Mesh,@MASS_Vol_diag_PWDG,omega);
% [I_mass_inn,J_mass_inn,M_inn] = assemMat_Inn_PDG2(Mesh,@MASS_Inn_PWDG,omega);
% [I_mass_bnd,J_mass_bnd,M_bnd] = assemMat_Bnd_PDG2(Mesh,[],@MASS_Bnd_PWDG,omega);
% I_mass = [I_mass_diag;I_mass_inn;I_mass_bnd];
% J_mass = [J_mass_diag;J_mass_inn;J_mass_bnd];
% M_mass = [M_diag;M_inn;M_bnd];
% M = sparse(I_mass,J_mass,M_mass);

% Assemble energy-norm inner product
[I_h1_diag,J_h1_diag,H_diag] = assemMat_Vol_PDG2(Mesh,@SCAPRO_Helm_diag_PWDG,omega);
[I_h1_inn,J_h1_inn,H_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega);
[I_h1_bnd,J_h1_bnd,H_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega);
I_h1 = [I_h1_diag;I_h1_inn;I_h1_bnd];
J_h1 = [J_h1_diag;J_h1_inn;J_h1_bnd];
H_h1 = [H_diag;H_inn;H_bnd];
B = sparse(I_h1,J_h1,H_h1);

% % Compute energy-norm projection of exact solution onto PWDG-space
% b_u = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega,u_ex,du_ex);
% u = B\b_u;

% Calculate error
err_p = EnergyErr_PWDG(Mesh,u,qr2,omega,u_ex,du_ex,0);
fprintf('Energy Error: %g    l2 norm: %g\n',err_p/omega,norm(u));

% Print coefficients
disp(abs(u));

% % Compute scalar products of pairs
% scapro = zeros(size(u));
% for j=1:ndir
%   ind = mod([j-1,j],ndir)+1;
%   scapro(j) = u(ind)'*B(ind,ind)*u(ind);
% end
% disp(abs(scapro));

% Compute projection of solution onto individual basis functions
proj1 = zeros(size(u));
for j=1:ndir
  w = B(j,j)\B(j,:)*u;
  proj1(j) = w'*B(j,j)*w;
end
% disp(abs(proj1));
[dummy,ind1] = sort(abs(proj1),1,'descend');

% % Iteratively project errors
% dec1 = zeros(size(u));
% u1 = zeros(size(u));
% for j=1:ndir
%   k = ind1(j);
%   w = B(k,k)\B(k,:)*(u-u1);
%   dec1(k) = w'*B(k,k)*w;
%   u1(k) = w;
% end
% err1 = sqrt(abs((u-u1)'*B*(u-u1)));
% fprintf('Error of Decomposition 1: %g    l2 norm: %g\n',err1/omega,norm(u1));
% disp(abs(u1));

% Iteratively solve local problem
u1a = zeros(size(u));
for j=1:ndir
  k = ind1(j);
  u1a(k) = A(k,k)\A(k,:)*(u-u1a);
end
err1a = sqrt(abs((u-u1a)'*B*(u-u1a)));
fprintf('Error of Decomposition 1a: %g    l2 norm: %g\n',err1a/omega,norm(u1a));
disp(abs(u1a));

% % Iteratively projected errors 2
% dec2 = zeros(size(u));
% u2 = zeros(size(u));
% ind2 = 1:ndir;
% for j=1:ndir
%   w = zeros(size(u));
%   for k=ind2
%     w(k) = B(k,k)\B(k,:)*(u-u2);
%   end
%   [w0,k0] = max(abs(w));
%   ind2 = setdiff(ind2,k0);
%   dec2(k0) = w(k0)'*B(k0,k0)*w(k0);
%   u2(k0) = w(k0);
% end
% err2 = sqrt(abs((u-u2)'*B*(u-u2)));
% fprintf('Error of Decomposition 2: %g    l2 norm: %g\n',err2/omega,norm(u2));
% disp(abs(u2));

% Iteratively solve local problem 2
u2a = zeros(size(u));
ind2a = 1:ndir;
for j=1:ndir
  w = zeros(size(u));
  for k=ind2a
    w(k) = A(k,k)\A(k,:)*(u-u2a);
  end
  [w0,k0] = max(abs(w));
  ind2a = setdiff(ind2a,k0);
  u2a(k0) = w(k0);
end
err2a = sqrt(abs((u-u2a)'*B*(u-u2a)));
fprintf('Error of Decomposition 2a: %g    l2 norm: %g\n',err2a/omega,norm(u2a));
disp(abs(u2a));

% % Iteratively projected errors 3
% u3 = zeros(size(u));
% n = ndir;
% for iter3=0:100
%   for j=1:n
%     w = zeros(size(u));
%     for k=1:ndir
%       w(k) = B(k,k)\B(k,:)*(u-u3);
%     end
%     [w0,k0] = max(abs(w));
%     u3(k0) = u3(k0) + w(k0);
%   end
%   err3 = sqrt(abs((u-u3)'*B*(u-u3)));
%   if(err3<err_p)
%     break;
%   end
% end
% fprintf('Error of Decomposition 3: %g    l2 norm: %g   (%g iterations)\n',err3/omega,norm(u3),iter3*n);
% disp(abs(u3));

% % Compute projection of solution onto pairs of basis functions
% proj2 = zeros(size(u));
% for j=1:ndir
%   ind = mod([j-1,j],ndir)+1;
%   w = B(ind,ind)\B(ind,:)*u;
%   proj2(j) = w'*B(ind,ind)*w;
% end
% disp(abs(proj2));

% % Solve regularized equation
% t = 1e-12;
% ut = (A'*A+t*speye(size(A,1)))\(A'*b);
% errt = sqrt(abs((u-ut)'*B*(u-u2)));
% fprintf('Error of Tychonov: %g    l2 norm: %g\n',errt/omega,norm(ut));
% disp(abs(ut));

return