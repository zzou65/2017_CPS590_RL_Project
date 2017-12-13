function [k,sv,v] = dispersion(ndir,omega,theta,flux_params,nVert)
%DISPERSION dispersion analysis of PWDG
%
%   K = DISPERSION(NDIR,OMEGA,THETA,FLUX_PARAMS,NVERT) attempts to
%   calculate K such that there exists a discrete solution of the Helmholtz
%   equation whose coefficients are periodic with wave number
%       K*(cos(THETA),sin(THETA)) .
%
%   NDIR is the number of local plane wave basis functions.
%
%   OMEGA is the wave number of the Helmholtz equation and of the plane
%   wave basis functions.
%
%   THETA is the angle to the x-axis of the wave vector of the plane wave.
%   It may be a vector.  If THETA is empty or not given, the maximal value
%   of k, sampled on a fine grid, is returned.
%
%   FLUX_PARAMS is a cell array containing the flux parameter names and
%   values.  The default values are equivalent to
%       FLUX_PARAMS = {'a',0.5,'b',0.5,'c',[0 0],'d',0.5} .
%
%   NVERT is the number of vertices in an element of the mesh, ie. 3 for
%   a triangular mesh (default) and 4 for a square mesh.
%
%   This code attempts to find the dispersion relation by minimizing the
%   smallest singular value of the stiffness matrix on an arbitrary element
%   of an infinite grid for periodic solutions.  This singular value is
%   returned by [K,SV] = DISPERSION(...).

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% % Define parameters for searching for k, ie kmin*omega < k < kmax*omega
% kmax = 2;
% kmin = 0.5;

% Read arguments
if(nargin<3 || isempty(theta))
  getmax = true;
  ntheta = ndir;
  theta = pi*(1/ndir:2/ndir:2-1/ndir);
else
  getmax = false;
  ntheta = numel(theta);
end
if(nargin<4 || isempty(flux_params))
  flux_params = {};
end
if(nargin<5 || isempty(nVert))
  nVert = 3;
end

if(nargout>=3)
  getev = true;
else
  getev = false;
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
% Mesh = orient_Elems(Mesh);
Mesh = add_Edge2Elem(Mesh);
Mesh = add_DGData(Mesh);
Mesh = set_Data_PWDG(Mesh,'Dir',dir,'Omega',omega,flux_params{:});
% plot_Mesh(Mesh,'petas');

% Determine indices for basis functions of each element
%   ind(j,:) are the indices for element j
ind = zeros(ndir,nElem);
ind(:) = 1:nElem*ndir;
ind = ind.';

% Determine position of element with respect to central elements
% trans(j,:) is the (row) vector to the basis point of the element j
trans = Mesh.Coordinates(Mesh.Elements(:,1),:);

% Determine indices of cetral elements corresponding to global indices
if(nVert==3)
  ind0 = ind([1 2 2 2 1 1],:);
  nCent = 2;
elseif(nVert==4)
  ind0 = ind(ones(1,nElem),:);
  nCent = 1;
end

% Initialize eigenvectors
if(getev)
  v = zeros(ndir*nCent,ntheta);
end

% Assemble stiffness matrix
[I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
A0 = accumarray([I_Inn,J_Inn],A_Inn);
A1 = A0(ind(1:nCent,:)',:);

% Construct stiffness matrix for k-periodic solution of central elements
A = @(k,theta) getA(k,theta,A1,ind0,ind,trans,nElem,nCent);

% Get smallest singular value of A
minS = @(k,theta) arrayfun(@(k1,k2)min(svd(A([k1,k2],theta))),k(:,1),k(:,2));

% % Define mesh of angles
% ntheta = 100;
% theta = 2*pi*(0:1/ntheta:1-1/ntheta);

% Find wave number for each angle
k = nan(size(theta));
sv = nan(size(theta));
% omega = complex(omega);
for j=1:ntheta
  [k0,sv(j)] = fminsearch(@(x) minS(x,theta(j)),[omega,0],...
    optimset('MaxFunEvals',2500,'MaxIter',1600,'TolX',1e-12));
  k(j) = complex(k0(1),k0(2));
%   [k(j),sv(j)] = fminsearch(@(x) minS(x,theta(j)),omega,optimset('TolX',1e-12));
%   [k(j),sv(j)] = fminbnd(@(x) minS(x,theta(j)),kmin*omega,kmax*omega,optimset('TolX',1e-12));

  if(getev) % compute eigenvector
    A_ = A(k(j),theta(j));
    [V,D] = eig(A_);
    d = diag(D);
    [ew,ew_ind] = min(abs(d));
    v(:,j) = V(:,ew_ind);
  end
end

% Get maximum if required
if(getmax)
  [dk,j] = max(abs(k-omega));
  k = k(j);
  sv = sv(j);
end

return

% construct periodic stiffness matrix depending on K
function A = getA(k,theta,A1,ind0,ind,trans,nElem,nCent)

if(length(k)>1)
  k = complex(k(1),k(2));
end

expK = @(j,k,theta) exp(i*trans(j,:)*k*[cos(theta);sin(theta)]);
A = A1(:,ind(1:nCent,:)');
for j=nCent+1:nElem
  A(:,ind0(j,:)) = A(:,ind0(j,:)) + expK(j,k,theta)*A1(:,ind(j,:));
end