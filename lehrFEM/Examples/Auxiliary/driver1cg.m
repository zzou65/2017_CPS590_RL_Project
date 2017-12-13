% Numerical experiment I: driver routine

clear all;

% Scaling parameter
tauvals = 10.^(-4:4);

% Load mesh
load 'PolygonMesh.mat';

% plot_Mesh(PolygonMesh,'a');

% Number of refinement steps
NREFS = 5;

itnum = zeros(NREFS+1,length(tauvals));

for k=1:length(tauvals)
  PolygonMesh.cg = true;
  fprintf('>>>> tau = %f <<<<<\n',tauvals(k));
  extev = auxspcexp1(1,tauvals(k),NREFS,PolygonMesh);
  itnum(:,k) = extev(:,1);
end

save 'PolygonMeshCGGS' itnum tauvals extev NREFS 
%save 'PolygonMeshCG' itnum tauvals extev NREFS 

clear all;

% Scaling parameter
tauvals = 10.^(-4:4);

% Load mesh
load 'Lshape.mat';

% plot_Mesh(Lshape,'a');

% Number of refinement steps
NREFS = 5;

cnd = zeros(NREFS+1,length(tauvals));

for k=1:length(tauvals)
  Lshape.cg = true;
  fprintf('>>>> tau = %f <<<<<\n',tauvals(k));
  extev = auxspcexp1(1,tauvals(k),NREFS,Lshape);
  cnd(:,k) = extev(:,1);
end

save 'LshapeCGGS' cnd tauvals extev NREFS 
% save 'LshapeCG' cnd tauvals extev NREFS 

