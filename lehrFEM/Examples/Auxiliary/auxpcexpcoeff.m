function extev = auxspcexpcoeff(Mesh,alpha,beta,NREFS)
% Numerical experiment studying nodal auxiliary space
% preconditioner for H(curl)-elliptic boundary value
% problems in the case of discontinuous coefficients
%
% Mesh -> mesh data structure with 'ElemFlag' field set
%
% alpha -> coefficient for curl-curl part for element flag ~= 1
% beta  -> coefficient for zero order part for element flag ~= 1
% (coefficients = 1, where element flag == 1)
% 
% NREFS = number of regular refinement steps
%

if (nargin < 4), NREFS = 5; end
if (nargin < 3), error('Insufficient arguments'); end
if (~isfield(Mesh,'ElemFlag')), error('Element flags missing'); end

% Switch: condition number computation <-> PCG
if (isfield(Mesh,'cg'))
  cgmon = true;
else
  cgmon = false;
end

ALPHA = @(x,elemflag,varargin) COEFF(x,elemflag,alpha);
BETA  = @(x,elemflag,varargin) COEFF(x,elemflag,beta);

fprintf('auxspcexpcoeff: Coefficient alpha in Omega_1 = %f\n',alpha); 
fprintf('auxspcexpcoeff: Coefficient beta in Omega_1 = %f\n',beta); 

% Prepare mesh for refinement
Mesh = add_Edges(Mesh);
Mesh = add_Edge2Elem(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

extev = [];
for ri=0:NREFS
  Mesh.verbose = true;
  
  fprintf('auxspcexpcoeff: Refinement level %d:\n',ri);
  if (cgmon)
    [itnum,flag,relres] = auxpccg(Mesh,ALPHA,BETA);
    fprintf('auxspcexpcoeff: itnum = %d, flag = %d\n',itnum,flag);
    extev = [extev; itnum flag relres size(Mesh.Elements,1) get_MeshMin(Mesh)];
  else
    [lmin,lmax,flag,iter] = auxpccond(Mesh,ALPHA,BETA);
    fprintf('auxspcexpcoeff: lmin = %f, lmax = %f [flag = %d, iter = %d]\n',...
	    lmin,lmax,flag,iter);
    extev = [extev; lmin lmax size(Mesh.Elements,1) get_MeshMin(Mesh)];
  end 
  Mesh = refine_REG(Mesh); 
end

end

% ======================================================================
% Local function providing coefficients 

function f = COEFF(x,elemflag,value)
% Coefficient for curl-curl-part of operator

f = ones(size(x,1),1);
if (elemflag ~= 1), f = value*f; end

end

