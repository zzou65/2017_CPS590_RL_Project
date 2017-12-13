function [A,P,VL,Grad,Lpot,f_load] = auxpcops(Mesh,ALPHA,BETA,F_HANDLE);
% Compute matrices required for nodal auxiliary space 
% preconditioning in H(curl)
%
% Mesh -> triangular mesh
%
% ALPHA,BETA -> coefficient functions
% F_HANDLE -> Source vector field @(x,elemflag,varargin) 
% (given a n x 2 matrix whose rows correspond to point coordinates,
%  this function must return another n x 2 matrix, whose rows give
% the values of the vectorfield at the points.
%
% See documentation below for output values
%

if (isa(BETA,'numeric'))
  beta = BETA(1); 
  BETA = @(x,varargin) beta*ones(size(x,1),1);
end

if (isa(ALPHA,'numeric'))
  alpha = ALPHA(1); 
  ALPHA = @(x,varargin) alpha*ones(size(x,1),1);
end

if (~isa(BETA,'function_handle')), error('BETA ~= numeric/function handle'); end
if (~isa(ALPHA,'function_handle')), error('ALPHA ~= numeric/fnct. handle'); end

if (~isfield(Mesh,'Edges'))
  Mesh = add_Edges(Mesh);
end

if (~isfield(Mesh,'Edge2Elem'))
  Mesh = add_Edge2Elem(Mesh);
end

% Extract boundary from mesh

BdEdges = get_BdEdges(Mesh);         %column vector of boundary edge indices
BdNodes = unique(Mesh.Edges(BdEdges,:)); % indices of boundary nodes
nNodes = size(Mesh.Coordinates,1);      
nEdges = size(Mesh.Edges,1);
IntNodes = setdiff(1:nNodes,BdNodes); % row vector
IntEdges = setdiff(1:nEdges,BdEdges); % row vector

if (isfield(Mesh,'verbose'))
  disp('Auxiliary space preconditioning for egde elements');
  disp('(Zero boundary conditions for auxiliary space)');
  fprintf('Mesh: %d nodes and %d edges\n',nNodes,nEdges);
  fprintf('Mesh: %d internal nodes, %d internal edges\n',...
	  length(IntNodes),length(IntEdges));
  fprintf('ALPHA'); disp(ALPHA);
  fprintf('BETA'); disp(BETA);
%   if (isa(beta,'numeric')), fprintf('beta = %f\n',beta); end
%   if (isa(alpha,'numeric')), fprintf('alpha = %f\n',alpha); end
end

% ======================================================================
% Assemble the stiffness matrix for edge element Galerkin
% discretization

% curl*alpha*curl-part of operator
[ICe,JCe,Ce] = assemMat_W1F(Mesh,@STIMA_Curl_W1F,ALPHA,P7O6());
% Assemble mass matrix for edge elements
[IMe,JMe,Me] = assemMat_W1F(Mesh,@MASS_W1F,BETA,P7O6());
% Build complete edge element Galerkin matrix A
A = sparse([ICe;IMe],[JCe;JMe],[Ce;Me]);
% Eliminate d.o.f. on the boundary
A = A(IntEdges,IntEdges); 

% ======================================================================
% Assemble the discrete Laplacian matrix used for correction in discrete
% potential space. This relies on the coefficient functions for the 
% zero order part of the operator

Lpot = assemMat_LFE(Mesh,@STIMA_Heat_LFE,P7O6(),BETA);
Lpot = Lpot(IntNodes,IntNodes);

% ======================================================================
% Assemble vector Laplacian/mass matrix for space of piecewise linear
% continuous vectorfields with zero boundary values.
% for 2nd order part use coefficient ALPHA
% for zero order part use coefficient BETA
% 
% NUMBERING CONVENTION (for d.o.f.s for "nodal' vectorfields
%  x-components: 1:nNodes
%  y-components: nNodes+1:2*nNodes

IntVNodes = [IntNodes,IntNodes+nNodes];  % row vector !

% Assemble block diagonal matrix
[IVL,JVL,VLe] = assemMat_LFE(Mesh,@STIMA_Heat_LFE,P7O6(),ALPHA);
[IVM,JVM,VMe] = assemMat_LFE(Mesh,@MASS_Weight_LFE,P7O6(),BETA);

VL = sparse([IVL;IVL+nNodes;IVM;IVM+nNodes],...
	    [JVL;JVL+nNodes;JVM;JVM+nNodes],...
	    [VLe;VLe;VMe;VMe]);
VL = VL(IntVNodes,IntVNodes);

% ======================================================================
% Assemble transfer matrices

% 1. discrete gradient matrix (edge-vertex incidence matrix)
Grad = Mat_G2E(Mesh); 
Grad = Grad(IntEdges,IntNodes);

% 2. edge interpolation matrix 
%    (nodal representation -> edge element representation)
P = convert_LFE2_W1F(Mesh);
P = P(IntEdges,IntVNodes);

% ======================================================================
% [Optional] assemble right hand side from source function
%

if ((nargin > 3) && (nargout > 5))
  f_load = assemLoad_W1F(Mesh,P7O6(),F_HANDLE);
  f_load = f_load(IntEdges);
else
  f_load = [];
end
