function [] = convstart_proj()
%CONVSTART_PROJ start of h-convergence

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Define problem parameters
omega = 2.^(3:7);
numOmega = numel(omega);
d = [1 0.5];
d = d/norm(d);
psi = 0;
lambda = 0.5*omega;
u_ex = @(x,k,varargin) sin(lambda(k)*x*d'+psi);
du_ex = @(x,k,varargin) lambda(k)*cos(lambda(k)*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);

% Define discretization parameters
nref = [1 5];
numRef = nref(2)-nref(1)+1;
ndir = 5;
phi = 2*pi*(0:1/ndir:1-1/ndir)';
dir = [cos(phi),sin(phi)];

% Define quadrature rules
% qr1 = gauleg(0,1,12);
qr2 = Duffy(TProd(gauleg(0,1,12))); % for triangular elements
% qr2 = TProd(gauleg(0,1,12)); % for quadrilateral elements

% Calculate mesh width
h = 2.^(1-nref(1):-1:1-nref(2));

% Initialize errors
err = nan(numOmega,numRef);

% Initialize mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
Mesh.Elements = [1 2 3; 1 4 3]; % triangular elements
% Mesh.Elements = [1 2 3 4]; % quadrilateral elements

Mesh = add_Edges(Mesh);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;
for j=1:nref(1)-1
  Mesh = refine_REG(Mesh);
end

% Loop over meshes
for j=1:numRef
  
  % Refine mesh and add fields
  Mesh = refine_REG(Mesh);
  Mesh = orient_Elems(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  Mesh = set_Data_PWDG(Mesh,'Dir',dir);
  
  % Loop over omegas
  for k=1:numOmega
    
    % Set omega
    Mesh = set_Data_PWDG(Mesh,'Omega',omega(k));
    
    % Assemble energy-norm inner product
    [I_h1_inn,J_h1_inn,H_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega(k));
    [I_h1_bnd,J_h1_bnd,H_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega(k));
    I_h1 = [I_h1_inn;I_h1_bnd];
    J_h1 = [J_h1_inn;J_h1_bnd];
    H_h1 = [H_inn;H_bnd];
    B = sparse(I_h1,J_h1,H_h1);

    % Calculate projection of exact solution onto discretization space
    b_proj = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega(k),u_ex,du_ex,k);
    u_proj = B\b_proj;
  
    % Calculate error
    err(k,j) = EnergyErr_PWDG(Mesh,u_proj,qr2,omega(k),u_ex,du_ex,0,k);
    
  end
  
end

% Plot errors
figure;
lgd = cell(size(omega));
for k=1:numOmega
  plot(h*omega(k),err(k,:));
  hold all;
  lgd{k} = sprintf('%s = %g, %s = %g','\omega',omega(k),'\lambda',lambda(k));
end
hold off;
set(gca,'XScale','log','XDir','reverse','YScale','log');
grid on;
legend(lgd{:},'Location','SouthWest');
xlabel('\bf h\omega');
ylabel('\bf energy-norm discretization error');
title(sprintf('%s PWDG discretization error for u(x) = sin(%sx%sd)','\bf','\lambda','\cdot'));