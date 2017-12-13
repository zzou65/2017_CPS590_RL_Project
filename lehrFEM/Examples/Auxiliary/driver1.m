% Numerical experiment I: driver routine

clear all;

% Scaling parameter
tauvals = 10.^(-4:4);

% Load mesh
load 'PolygonMesh.mat';

% plot_Mesh(PolygonMesh,'a');

% Number of refinement steps
NREFS = 5;

cnd = zeros(NREFS+1,length(tauvals));

for k=1:length(tauvals)
  fprintf('>>>> tau = %f <<<<<\n',tauvals(k));
  extev = auxspcexp1(1,tauvals(k),NREFS,PolygonMesh);
  cnd(:,k) = extev(:,2)./extev(:,1);
end

save 'PolygonMeshCndGS' cnd tauvals extev NREFS 
%save 'PolygonMeshCnd' cnd tauvals extev NREFS 

% figure('name','PolygonMesh');
% [X,Y] = meshgrid(0:NREFS,tauvals);
% mesh(X,Y,cnd);
% title('{\bf Convex polygon}');
% set(gca,'xscale','log');
% xlabel('{\bf \tau}','fontsize',14);
% ylabel('{\bf refinement level}');
% zlabel('{\bf spectral condition number}');

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
  fprintf('>>>> tau = %f <<<<<\n',tauvals(k));
  extev = auxspcexp1(1,tauvals(k),NREFS,Lshape);
  cnd(:,k) = extev(:,2)./extev(:,1);
end

save 'LshapeCndGS' cnd tauvals extev NREFS 
% save 'LshapeCnd' cnd tauvals extev NREFS 

% figure('name','Lshape');
% [X,Y] = meshgrid(0:NREFS,tauvals);
% mesh(X,Y,cnd);
% title('{\bf L-shaped domain}');
% set(gca,'xscale','log');
% xlabel('{\bf \tau}','fontsize',14);
% ylabel('{\bf refinement level}');
% zlabel('{\bf spectral condition number}');
