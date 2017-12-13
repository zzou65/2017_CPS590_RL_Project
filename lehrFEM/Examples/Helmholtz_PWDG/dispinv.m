function omega = dispinv(ndir,k,flux_params,nVert,nOmega)
%DISPINV disperion analysis of PWDG
%   
%   OMEGA = DISPINV(NDIR,K,PARAMS,NVERT,NOMEGA) calculates, for each
%   row k of K, a value OMEGA such that the discretization using plane 
%   waves with wave number |k| and periodic with wave vector k satisfies
%   the Helmholtz equation with wave number OMEGA.
%
%   NDIR is the number of local plane wave basis functions.
%
%   K is a N-by-2 matrix; OMEGA is calculated for each row of K.
%
%   FLUX_PARAMS is a cell array containing the flux parameter names and
%   values.  The default values are equivalent to
%       FLUX_PARAMS = {'a',0.5,'b',0.5,'c',[0 0],'d',0.5} .
%
%   NVERT is the number of vertices in an element of the mesh, ie. 3 for
%   a triangular mesh (default) and 4 for a square mesh.
%
%   NOMEGA is the number of values for OMEGA to return per row of K.  The
%   default value is 1.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Read arguments
if(nargin<3 || isempty(flux_params))
  flux_params = {};
end
if(nargin<4 || isempty(nVert))
  nVert = 3;
end
if(nargin<5 || isempty(nOmega))
  nOmega = 1;
end

% Define discretization parameters
phi = 2*pi*(0:1/ndir:1-1/ndir)';
dir = [cos(phi),sin(phi)];

% Generate mesh
if(nVert==3)
  Mesh.Coordinates = [0 0; 1 0; 1 1; 0 1; 0 -1; 2 1; 1 2; -1 0];
  Mesh.Elements = [1 2 3; 1 3 4; 5 2 1; 2 6 3; 4 3 7; 8 1 4];
elseif(nVert==4)
  Mesh.Coordinates = [0 0; 1 0; 1 1; 0 1;...
    0 -1; 1 -1; 2 0; 2 1; 1 2; 0 2; -1 1; -1 0];
  Mesh.Elements = [1 2 3 4; 5 6 2 1; 2 7 8 3; 4 3 9 10; 12 1 4 11];
end
nElem = size(Mesh.Elements,1);
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;
Mesh = add_Edge2Elem(Mesh);
Mesh = add_DGData(Mesh);
Mesh = set_Data_PWDG(Mesh,'Dir',dir);

% Determine indices for basis functions of each element
%   ind(j,:) are the indices for element j
ind = zeros(ndir,nElem);
ind(:) = 1:nElem*ndir;
ind = ind.';

% Determine position of element with respect to central elements
%   trans(j,:) is the (row) vector to the basis point of the element j
trans = Mesh.Coordinates(Mesh.Elements(:,1),:);

% Determine indices of cetral elements corresponding to global indices
if(nVert==3)
  ind0 = ind([1 2 2 2 1 1],:);
  nCent = 2;
elseif(nVert==4)
  ind0 = ind(ones(1,nElem),:);
  nCent = 1;
end

% Correct value of nOmega if out of range
nOmega = max(1,min(nCent*ndir,nOmega));

% Initialize solution omega
nk = size(k,1);
omega = nan(nk,nOmega);

% Loop over values of k
for j=1:nk
  
  % Set flux parameters
  k_ = norm(k(j,:));
  Mesh = set_Data_PWDG(Mesh,'Omega',k_,flux_params{:});
  
  % Assemble stiffness matrix
  [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,k_);
  A0 = accumarray([I_Inn,J_Inn],A_Inn);
  A1 = A0(ind(1:nCent,:)',:);
  
  % Assemble mass matrix
  [I_mass_inn,J_mass_inn,M_inn] = assemMat_Inn_PDG2(Mesh,@MASS_Inn_PWDG,k_);
  M0 = accumarray([I_mass_inn,J_mass_inn],M_inn);
  M1 = M0(ind(1:nCent,:)',:);
  
  % Construct stiffness and mass matrices for k-periodic solution of central elements
  expK = exp(i*trans*k(j,:)');
  A = A1(:,ind(1:nCent,:)');
  M = M1(:,ind(1:nCent,:)');
  for l=nCent+1:nElem
    A(:,ind0(l,:)) = A(:,ind0(l,:)) + expK(l)*A1(:,ind(l,:));
    M(:,ind0(l,:)) = M(:,ind0(l,:)) + expK(l)*M1(:,ind(l,:));
  end
  
  % Construct stiffness matrix for Laplacian
  A_ = A + k_^2*M;
  
  % Solve eigenvalue equation
  ev = eig(A_,M);
%   [d,omega_ind] = sort(abs(imag(sqrt(ev)-k_))); % is 'imag' reasonable?
%   [d,omega_ind] = sort(abs(real(sqrt(ev)-k_))); % very noisy...
%   [d,omega_ind] = min(abs(sqrt(ev)-k_));
  [dk,omega_ind] = sort(abs(sqrt(ev)-k_));
  omega(j,:) = sqrt(ev(omega_ind(1:nOmega)));

end