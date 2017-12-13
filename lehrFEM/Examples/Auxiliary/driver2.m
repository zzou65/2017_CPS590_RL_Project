% Numerical experiment II: driver routine

clear all;

% Scaling parameter
tauvals = 10.^[-4,-2,0,2,4];

% Load mesh
load 'Lshape.mat';

% Lshape.plot = [5 10 20 40];

% Number of refinement steps
NREFS = 35;

cnd = [];

for k=1:length(tauvals)
  fprintf('>>>> tau = %f <<<<<\n',tauvals(k));
  extev = auxspcexp1(1,tauvals(k),NREFS,Lshape,[-0.1,0.1,-0.1,0.1]);
  cnd = [cnd, extev(:,2)./extev(:,1)];
end

save 'LshapeCndlocal' cnd tauvals extev NREFS 

% figure('name','Lshape');
% [X,Y] = meshgrid(0:NREFS,tauvals);
% mesh(X,Y,cnd);
% title('{\bf L-shaped domain}');
% set(gca,'xscale','log');
% xlabel('{\bf \tau}','fontsize',14);
% ylabel('{\bf refinement level}');
% zlabel('{\bf spectral condition number}');
