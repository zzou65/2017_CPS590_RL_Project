function extev = auxspcexp1(alpha,beta,NREFS,Mesh,BBox)
% Numerical experiment: nodal auxiliary space
% preconditioning for edge elements
% on a sequence of regularly refined meshes of the
% unit square
%
% alpha -> coefficient for curl-curl-part
% beta  -> coefficient for zero order part
% NREFS -> Number of refinement steps
% BBox -> bounding box controlling local refinement
% (only cells with a vertex inside BBox are refined)
%
% RETURN VALUE:
%
% NREFSx4-matrix:
% column 1 -> minimal eigenvalue
% column 2 -> maximal eigenvalue
% column 3 -> #elements 
% column 4 -> minimal edge length

if (nargin < 3), NREFS = 5; end
if (nargin < 2), beta = 1; end
if (nargin < 1), alpha = 1; end;

MAXCELLS = 100000; % Maximum number of cells for adaptive refinement
BBSCHRINK = 0.8;  % Shrinking factor for bounding box

if (nargin < 4)
  Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
  Mesh.Elements = [1 2 4;2 3 4];
end

Mesh = add_Edges(Mesh);
Mesh = add_Edge2Elem(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

if (nargin < 4), Mesh = refine_REG(Mesh); end

% Field specifying on which level the mesh should be plotted
if (isfield(Mesh,'plot')) 
  plotlevs = Mesh.plot; 
else 
  plotlevs = []; 
end

% Switch: condition number computation <-> PCG
if (isfield(Mesh,'cg'))
  cgmon = true;
else
  cgmon = false;
end

extev = [];
for ri=0:NREFS
  Mesh.verbose = true;
  
  fprintf('Refinement level %d:\n',ri);
  if (cgmon)
    [itnum,flag,relres] = auxpccg(Mesh,alpha,beta);
    fprintf('itnum = %d, flag = %d\n',itnum,flag);
    extev = [extev; itnum flag relres size(Mesh.Elements,1) get_MeshMin(Mesh)];
  else
    [lmin,lmax,flag,iter] = auxpccond(Mesh,alpha,beta);
    fprintf('lmin = %f, lmax = %f [flag = %d, iter = %d]\n',...
	    lmin,lmax,flag,iter);
    extev = [extev; lmin lmax size(Mesh.Elements,1) get_MeshMin(Mesh)];
  end

  if (nargin < 5), Mesh = refine_REG(Mesh); 
  else
    fprintf('Bounding box = '); disp(BBox);
    Mesh = init_LEB(Mesh);
    marked = [];
    for i=1:size(Mesh.Elements,1)
      flag = false;
      for k=1:3
	P = Mesh.Coordinates(Mesh.Elements(i,k),:);
	if ((BBox(1) < P(1)) && ...
	    (BBox(2) > P(1)) && ...
	    (BBox(3) < P(2)) && ...
	    (BBox(4) > P(2))), flag = true;
	end
      end
      if (flag), marked = [marked,i]; end
    end
    Mesh = refine_LEB(Mesh,marked);
    Mesh = add_Edges(Mesh);
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1; 
    
    nCells = size(Mesh.Coordinates,1);
    fprintf('REFINEMENT STEP %d: %d cells\n',ri,nCells);

    if (sum(plotlevs == ri) > 0)
      Plot_Mesh(Mesh,'a'); 
      title(sprintf('{\\bf Locally refined mesh level %d: h_{min}=%g}',...
		    ri,get_MeshMin(Mesh)));
      drawnow;
      print('-depsc2',sprintf('reflev%d.eps',ri));
    end
    
    if (nCells > MAXCELLS), break; end
    BBox = 0.75*BBox;
  end
end




