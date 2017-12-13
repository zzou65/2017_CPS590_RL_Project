function [fig1,fig2,fig3] = pollut_omega(refs,ndir,nvert,omega,d,flux_params,mtd_name)
%POLLUT_OMEGA plot pollution error vs. omega
%   
%   POLUT_OMEGA(REFS,NDIR,NVERT,OMEGA,D,FLUX_PARAMS,MTD_NAME) plots
%   the pollution error of the discontinuous plane wave discretization of
%   the Helholtz equation with wave number OMEGA against OMEGA.  The
%   pollution error is the error in the discrete solution minus the
%   projection error, divided by the projection error and maximized over
%   the mesh width.  It measures the maximal discrepency between the
%   convergence of the discrete solution and the best approximation.
%
%   This code also plots estimates for start of convergence based on the
%   maximal pollution error.  A value for h*omega is determined for which
%   the projection error plus the maximal polution error is equal to the
%   total error on coarse meshes.
%
%   Additionaly, the convergence for all omega is plotted together with the
%   projection error and the estimates for the start of convergence.
%
%   [FIG1,FIG2,FIG3] = POLLUT_OMEGA(...) returns pointers to the three
%   figures.
%
%   REFS is a vector of length 2 specifying the number of refinements on
%   the coarsest and finest mesh.
%
%   NDIR is the number of local plane wave basis functions.
%
%   NVERT is the number of vertices in an element of the mesh, ie. 3 for
%   a triangular mesh (default) and 4 for a square mesh.
%
%   D is the propagation direction of the exact solution in the sample
%   problem.  If D is empty, a value with large dispersion is chosen.
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
  refs = [1 4];
end
if(nargin<2 || isempty(ndir))
  ndir = 4;
end
if(nargin<3 || isempty(nvert))
  nvert = 3;
end
if(nargin<4 || isempty(omega))
  omega = 2.^(0:4);
end
if(nargin<6 || isempty(flux_params))
  flux_params = {};
end
if(nargin<7 || isempty(mtd_name))
  mtd_name = 'PWDG';
end
if(nargin<5 || isempty(d))
  % Determine d with large dispersion
%   d = [cos(pi/ndir),sin(pi/ndir)];
  theta = pi*(1/ndir:2/ndir:2-1/ndir)';
  k = [cos(theta),sin(theta)];
  disp = abs(dispinv(ndir,k,flux_params,nvert)-1);
  [disp_max,ind_max] = max(disp);
  d = k(ind_max,:);
end
nomega = numel(omega);
nref = refs(2)-refs(1)+1;           % number of mesh refinements
h = 2.^(1-refs(1):-1:1-refs(2));    % mesh width

% define dispersive sample problem
phi = 2*pi*(0:1/ndir:1-1/ndir)';    
dir = [cos(phi),sin(phi)];          % directions of plane wave basis functions
d = d/norm(d);
lambda = omega;
u_ex = @(x,k,varargin) exp(i*lambda(k)*x*d');
du_ex = @(x,k,varargin) i*lambda(k)*exp(i*lambda(k)*x*d([1 1],:)').*d(ones(size(x,1),1),:);
g = @(x,n,k,varargin) i*(lambda(k)*(d*n')+omega(k))*exp(i*lambda(k)*x*d');
f = @(x,k,varargin) (lambda(k)^2-omega(k)^2)*exp(i*lambda(k)*x*d');

% initialize quadrature rules
qr1 = gauleg(0,1,12);
if(nvert==3) % triangular elements
  qr2 = Duffy(TProd(gauleg(0,1,12)));
elseif(nvert==4) % quadrilateral elements
  qr2 = TProd(gauleg(0,1,12));
end

% initialize errors
err_proj = nan(nomega,nref);
err_pwdg = nan(nomega,nref);
err_pol = nan(nomega,nref);
homega = nan(nomega,nref);

% initialize mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
if(nvert==3) % triangular elements
  Mesh.Elements = [1 2 3; 1 4 3];
elseif(nvert==4) % quadrilateral elements
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
  Mesh = set_Data_PWDG(Mesh,'Dir',dir);

  % loop over values of omega
  for k=1:nomega
    
    homega(k,j) = h(j)*omega(k);
    
    % set flux parameters
    Mesh = set_Data_PWDG(Mesh,'Omega',omega(k),flux_params{:});

    % assemble stiffness matrix
    [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega(k));
    [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega(k));
    A = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

    % assemble load vector
    b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega(k),g,k);
    b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega(k),f,k);
    b = b_f + b_g;

    % solve equation
%     u = A\b;
    [PL,PR] = assemPrec_SVD_PDG(Mesh,A);
    w = (PL*A*PR)\(PL*b);
    u = PR*w;
    
    % assemble energy-norm inner product
    [I_h1_inn,J_h1_inn,H_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega(k));
    [I_h1_bnd,J_h1_bnd,H_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega(k));
    I_h1 = [I_h1_inn;I_h1_bnd];
    J_h1 = [J_h1_inn;J_h1_bnd];
    H_h1 = [H_inn;H_bnd];
    B = sparse(I_h1,J_h1,H_h1);

    % calculate projection of exact solution onto discretization space
    b_proj = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega(k),u_ex,du_ex,k);
    u_proj = B\b_proj;

    % calculate polution error
    err_proj(k,j) = EnergyErr_PWDG(Mesh,u_proj,qr2,omega(k),u_ex,du_ex,0,k)/omega(k);
    err_pwdg(k,j) = EnergyErr_PWDG(Mesh,u,qr2,omega(k),u_ex,du_ex,0,k)/omega(k);
    err_pol(k,j) = abs(err_pwdg(k,j)-err_proj(k,j))/err_proj(k,j);
%     du = u-u_proj;
%     err_pol(k,j) = sqrt(abs(du'*B*du))/(omega(k)*err_proj(k,j));

  end

end

% determine maximal polution error for each omega
[err_max,ind_max] = max(err_pol,[],2);
homega_max = h(ind_max).*omega;
err_max_full = nan(size(err_max));
for k=1:nomega
  err_max_full(k) = err_pwdg(k,ind_max(k));
end

% construct linear bound for convergence through the point with maximal
% deviation from the projection error, with slope equal to the slope of the
% projection error and cut off at the maximal value of the projection
% error.  The value of h*omega at this cutoff point is a measure for where
% the convergence begins.
log_err_proj = log(err_proj);
log_homega = log(homega);
err_0 = exp(max(log_err_proj(:,1)));      % initial error
d_log_err = diff(log_err_proj(:,end-1:end),1,2);
d_log_homega = diff(log_homega(:,end-1:end),1,2);
m_log = max(d_log_err./d_log_homega);     % slope of projection error in log-log plot
homega_0 = homega_max.*(err_0./err_max_full.').^(1/m_log); % start of convergence

% plot polution effect
fig1 = figure;
plot(omega,err_max);
set(gca,'XScale','log','XLim',[min(omega),max(omega)],'XTick',omega);
grid on;
xlabel('\bf \omega');
% ylabel('\bf max_h || u_{PWDG} - u_{PROJ} ||_\omega / || u_{PROJ} - u ||_\omega');
ylabel('\bf max_h | || err_{PWDG} ||_\omega - || err_{PROJ} ||_\omega | / || err_{PROJ} ||_\omega');
title(sprintf('\\bf Polution effect in %s, p=%g',mtd_name,ndir));

% plot estimate for start of convergence
fig2 = figure;
plot(omega,homega_0);
set(gca,'XScale','log','XLim',[min(omega),max(omega)],'XTick',omega,'YScale','log');
grid on;
xlabel('\bf \omega');
ylabel('\bf h\omega at which convergence begins (estimated)');
title(sprintf('\\bf Start of convergence of %s, p=%g',mtd_name,ndir));

p = polyfit(log(omega),log(homega_0),1);
add_Slope(gca,'SouthWest',p(1));

% plot bound for convergence of PWDG solution
fig3 = figure;
lgd = cell(1,nomega);
handle1 = nan(nomega,1);
hold on;
for k=1:nomega
  handle1(k) = plot(h*omega(k),err_pwdg(k,:));
  hold all;
  lgd{k} = sprintf('\\omega = %g',omega(k));
end
hold on;
for k=1:nomega
  plot(h*omega(k),err_proj(k,:),':');
  hold all;
end
hold on;
set(gca,'XScale','log','XDir','reverse','YScale','log');
homega_lim = get(gca,'XLim');
for k=1:nomega
  homega = [homega_lim(1),homega_0(k),homega_lim(2)];
  err_est = [err_max_full(k)*(homega_lim(1)/homega_max(k))^m_log,err_0,err_0];
  plot(homega,err_est,'--');
  hold all;  
end
hold on;
for k=1:nomega
  plot(homega_0(k),err_0,'o');
  hold all;  
end
handle2 = plot(nan,nan,'k:',nan,nan,'k--',nan,nan,'ko');
hold off;
grid on;
legend([handle1;handle2],lgd{:},'projection error','estimate','start of convergence','Location','SouthWest');
xlabel('\bf h\omega');
ylabel('\bf ( \omega^{-2}| error |_{1,h}^2 + || error ||_0^2 )^{1/2}');
title(sprintf('\\bf Estimate for convergence of %s, p=%g',mtd_name,ndir));
  

return