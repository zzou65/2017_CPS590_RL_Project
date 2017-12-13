function [] = conv_h_mixed()
%CONV_H_MIXED h-convergence for mixed method

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Define problem parameters
omega = 1;
d = [1 0.5];
d = d/norm(d);
lambda = omega;
psi = 0;
u_ex = @(x,varargin) sin(lambda*x*d'+psi);
du_ex = @(x,varargin) lambda*cos(lambda*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);
% gI = @(x,n,varargin) lambda*(d*n')*cos(lambda*x*d'+psi) + i*omega*sin(lambda*x*d'+psi);
gD = @(x,varargin) sin(lambda*x*d'+psi);
f = @(x,varargin) (lambda^2-omega^2)*sin(lambda*x*d'+psi);

% Define discretization parameters
ndir = 5;
nref = [0 3];
flux_params = {'a',@(w,h,varargin) 2/(w*h),'b',@(w,h,varargin) 0.01*w*h};

% Define quadrature rules
qr1 = gauleg(0,1,12);
qr2 = Duffy(TProd(gauleg(0,1,12))); % for triangular elements
% qr2 = TProd(gauleg(0,1,12)); % for quadrilateral elements

% Calculate mesh width
h = 2.^(1-nref(1):-1:1-nref(2));

% Initialize errors
num = nref(2)-nref(1)+1;
err_primal = nan(1,num);
err_mixed = nan(1,num);
err_proj = nan(1,num);

% Initialize mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
Mesh.Elements = [1 2 3; 1 4 3]; % triangular elements
% Mesh.Elements = [1 2 3 4]; % quadrilateral elements

Mesh = add_Edges(Mesh);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;
for j=1:nref(1)
  Mesh = refine_REG(Mesh);
end

% Loop over meshes
for j=1:num
  
  % Refine mesh and add fields
  Mesh = orient_Elems(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  Mesh = set_Data_PWDG(Mesh,'nDofs',ndir,'Omega',omega,flux_params{:});

  % Assemble stiffness matrix for primal method
  [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
%   [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega);
  [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Dir_Bnd_PWDG,omega);
  A_primal = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

  % Assemble load vector for primal method
%   b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
  b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Dir_Bnd_PWDG,qr1,omega,gD);
  b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,f);
  b_primal = b_f + b_g;
  
  % Solve equation for primal method
  u_primal = A_primal\b_primal;
  
  % Assemble L2 inner product
  [I_mass_inn,J_mass_inn,M_inn] = assemMat_Inn_PDG2(Mesh,@MASS_Inn_PWDG,omega);
  [I_mass_bnd,J_mass_bnd,M_bnd] = assemMat_Bnd_PDG2(Mesh,[],@MASS_Bnd_PWDG,omega);
  I_mass = [I_mass_inn;I_mass_bnd];
  J_mass = [J_mass_inn;J_mass_bnd];
  M_mass = [M_inn;M_bnd];
  M = sparse(I_mass,J_mass,M_mass);
  
  % Assemble stiffness matrix for mixed method
  [I_Vol,J_Vol,A_Vol] = assemMat_Vol_PDG2_vec(Mesh,3,@STIMA_Helm_Vol_PWDG_dual,omega,M);
  [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2_vec(Mesh,3,@STIMA_Helm_Inn_PWDG_dual,omega);
%   [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2_vec(Mesh,[],3,@STIMA_Helm_Imp_Bnd_PWDG_dual,omega);
  [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2_vec(Mesh,[],3,@STIMA_Helm_Dir_Bnd_PWDG_dual,omega);
  A_mixed = sparse([I_Vol;I_Inn;I_Bnd],[J_Vol;J_Inn;J_Bnd],[A_Vol;A_Inn;A_Bnd]);
  
  % Assemble load vector for primal method
%   b_g = assemLoad_Bnd_PDG2_vec(Mesh,[],3,@LOAD_Imp_Bnd_PWDG_dual,qr1,omega,gI);
  b_g = assemLoad_Bnd_PDG2_vec(Mesh,[],3,@LOAD_Dir_Bnd_PWDG_dual,qr1,omega,gD);
  b_f = assemLoad_Vol_PDG2_vec(Mesh,3,@LOAD_Vol_PWDG_dual,qr2,omega,f);
  b_mixed = b_f + b_g;
  
  % Solve equation for primal method
  U_mixed = A_mixed\b_mixed;
  u_mixed = U_mixed(1:end/3);
  
  % Assemble energy-norm inner product
  [I_h1_inn,J_h1_inn,H_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega);
  [I_h1_bnd,J_h1_bnd,H_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega);
  I_h1 = [I_h1_inn;I_h1_bnd];
  J_h1 = [J_h1_inn;J_h1_bnd];
  H_h1 = [H_inn;H_bnd];
  B = sparse(I_h1,J_h1,H_h1);
  
  % Calculate projection of exact solution onto discretization space
  b_proj = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega,u_ex,du_ex);
  u_proj = B\b_proj;
  
  % Calculate errors
  err_primal(j) = EnergyErr_PWDG(Mesh,u_primal,qr2,omega,u_ex,du_ex,0);
  err_mixed(j) = EnergyErr_PWDG(Mesh,u_mixed,qr2,omega,u_ex,du_ex,0);
  err_proj(j) = EnergyErr_PWDG(Mesh,u_proj,qr2,omega,u_ex,du_ex,0);
  
  % Refine mesh
  if(j<num)
    Mesh = refine_REG(Mesh);
  end
    
end

% Check if exact solution in kernel
kernstring = '';
if(omega==lambda)
  kernstring = ' for ex. sol. in kernel';
end
  
% Plot errors
figure;
plot(h*omega,err_primal,'-',h*omega,err_mixed,'--',h*omega,err_proj,':');
set(gca,'XScale','log','XDir','reverse','YScale','log','XLim',omega*[min(h),max(h)]);
grid on;
xlabel('\bf h\omega');
ylabel('\bf energy-norm error');
title(sprintf('\\bf Convergence of PWDG%s, %d local basis fn.',kernstring,ndir));
legend('primal','mixed','proj.');

% Add convergence rate
p = polyfit(log(h),log(err_proj),1);
add_Slope(gca,'SouthEast',p(1));
