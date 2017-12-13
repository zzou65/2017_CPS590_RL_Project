function fig = pollut_h(refs,nVert,flux_params,mtd_name)
%POLLUT_H pollution effect
%   
%   POLUT_H(REFS,NVERT,FLUX_PARAMS,MTD_NAME) plots the mesh width h
%   and the wave number omega of the Helmholtz equation such that for each
%   h, omega maximizes the pollution error.  The pollution error is the 
%   energy norm of the discrete solution minus the best approximation,
%   divided by the projection error and maximized over the mesh width.  It
%   measures the maximal discrepency between the convergence of the
%   discrete solution and the best approximation.  The dependence of
%   h*omega on omega indicates the extra mesh refinements required to
%   compensate for dispersion.
%
%   REFS is a vector of length 2 specifying the number of refinements on
%   the coarsest and finest mesh.
%
%   NVERT is the number of vertices in an element of the mesh, ie. 3 for
%   a triangular mesh (default) and 4 for a square mesh.
%
%   FLUX_PARAMS is a cell array containing the flux parameter names and
%   values.  The default values are equivalent to
%       FLUX_PARAMS = {'a',0.5,'b',0.5,'c',[0 0],'d',0.5} .
%
%   MTD_NAME is a string containing the name of the method used, for
%   example 'PWDG'.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% read arguments
if(nargin<1 || isempty(refs))
  refs = [1 3];
end
if(nargin<2 || isempty(nVert))
  nVert = 3;
end
if(nargin<3 || isempty(flux_params))
  flux_params = {};
end
if(nargin<4 || isempty(mtd_name))
  mtd_name = 'PWDG';
end
nref = refs(2)-refs(1)+1;           % number of mesh refinements
h = 2.^(1-refs(1):-1:1-refs(2));    % mesh width

% initialize quadrature rules
qr1 = gauleg(0,1,12);
if(nVert==3) % triangular elements
  qr2 = Duffy(TProd(gauleg(0,1,12)));
elseif(nVert==4) % quadrilateral elements
  qr2 = TProd(gauleg(0,1,12));
end

% initialize omega
omega = nan(1,nref);

% initialize mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
if(nVert==3) % triangular elements
  Mesh.Elements = [1 2 3; 1 4 3];
elseif(nVert==4) % quadrilateral elements
  Mesh.Elements = [1 2 3 4];
end
Mesh = add_Edges(Mesh);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;
for j=1:refs(1)-1
  Mesh = refine_REG(Mesh);
end

% loop over meshes
for j=1:nref

  % refine mesh and add fields
  Mesh = refine_REG(Mesh);
  Mesh = orient_Elems(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);

  % define function to minimize
  err_h = @(omega) -getErr(omega,Mesh,flux_params,qr1,qr2);
  
%   % define initial guess
%   omega_0 = 1;
  
  % define bounds for omega
  omega_min = 0.1/h(j);
  omega_max = 10/h(j);

  % minimize
%   omega(j) = fminsearch(err_h,omega_0);
  omega(j) = fminbnd(err_h,omega_min,omega_max);

end

% plot h vs. omega
fig = figure;
plot(omega,h.*omega);
set(gca,'XScale','log','YScale','log','XLim',[min(omega),max(omega)],'YDir','reverse');
grid on;
xlabel('\bf \omega');% = argmax || u_{PWDG} - u_{PROJ} ||_\omega / || u_{PROJ} - u ||_\omega');
ylabel('\bf h\omega');
title(sprintf('\\bf Maximal polution error for %s',mtd_name));

return

%%% calculate error for given omega
function err = getErr(omega,Mesh,flux_params,qr1,qr2)

% define dispersive sample problem
ndir = 4;
phi = 2*pi*(0:1/ndir:1-1/ndir)';    
dir = [cos(phi),sin(phi)];          % directions of plane wave basis functions
d = [1 1];
d = d/norm(d);
lambda = omega;
psi = 0;
u_ex = @(x,varargin) sin(lambda*x*d'+psi);
du_ex = @(x,varargin) lambda*cos(lambda*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);
g = @(x,n,varargin) lambda*(d*n')*cos(lambda*x*d'+psi) + i*omega*sin(lambda*x*d'+psi);
f = @(x,varargin) (lambda^2-omega^2)*sin(lambda*x*d'+psi);

% set PWDG data
Mesh = set_Data_PWDG(Mesh,'Dir',dir,'NewDir','replace','Omega',omega,flux_params{:});

% assemble stiffness matrix
[I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
[I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega);
A = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

% assemble load vector
b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,g);
b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,f);
b = b_f + b_g;

% solve equation
u = A\b;

% assemble energy-norm inner product
[I_h1_inn,J_h1_inn,H_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega);
[I_h1_bnd,J_h1_bnd,H_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega);
I_h1 = [I_h1_inn;I_h1_bnd];
J_h1 = [J_h1_inn;J_h1_bnd];
H_h1 = [H_inn;H_bnd];
B = sparse(I_h1,J_h1,H_h1);

% calculate projection of exact solution onto discretization space
b_proj = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega,u_ex,du_ex);
u_proj = B\b_proj;

% calculate polution error
err_proj = EnergyErr_PWDG(Mesh,u_proj,qr2,omega,u_ex,du_ex,0);
du = u-u_proj;
err = sqrt(abs(du'*B*du))/err_proj;

return
