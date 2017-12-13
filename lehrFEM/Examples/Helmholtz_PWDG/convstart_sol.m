function [] = convstart_sol()
%CONVSTART_SOL start of h-convergence

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
% lambda = ones(size(omega));
u_ex = @(x,k,varargin) sin(lambda(k)*x*d'+psi);
du_ex = @(x,k,varargin) lambda(k)*cos(lambda(k)*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);
g = @(x,n,k,varargin) lambda(k)*(d*n')*cos(lambda(k)*x*d'+psi) + i*omega(k)*sin(lambda(k)*x*d'+psi);
f = @(x,k,varargin) (lambda(k)^2-omega(k)^2)*sin(lambda(k)*x*d'+psi);

% Define discretization parameters
nref = [1 5];
numRef = nref(2)-nref(1)+1;
ndir = 5;
phi = 2*pi*(0:1/ndir:1-1/ndir)';
dir = [cos(phi),sin(phi)];

% Define quadrature rules
qr1 = gauleg(0,1,12);
qr2 = Duffy(TProd(gauleg(0,1,12))); % for triangular elements
% qr2 = TProd(gauleg(0,1,12)); % for quadrilateral elements

% Define flux parameters
flux_params = {'a',@(omega,h,varargin) 1/(h*omega),...
  'b',@(omega,h,varargin) 0.01*h*omega,'d',0.5};
% flux_params = {};  % use ultra-weak variational formulation

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
    
    % Set PWDG flux parameters
    Mesh = set_Data_PWDG(Mesh,'Omega',omega(k),flux_params{:});

    % Assemble stiffness matrix
    [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega(k));
    [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega(k));
    A = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

    % Assemble load vector
    b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega(k),g,k);
    b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega(k),f,k);
    b = b_f + b_g;

    % Solve equation
    u = A\b;
  
    % Calculate error
    err(k,j) = EnergyErr_PWDG(Mesh,u,qr2,omega(k),u_ex,du_ex,0,k);
    
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
ylabel('\bf energy-norm error');
title(sprintf('%s Convergence of PWDG for u(x) = sin(%sx%sd)','\bf','\lambda','\cdot'));