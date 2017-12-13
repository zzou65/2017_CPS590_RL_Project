% Numerical experiment I: driver routine

clear all;

% Scaling parameter
coeffvals = [0.001,0.01,0.1,2,5,10,20,50,100,200,1000];

% Load mesh
load 'TriangleSubDom.mat';

% Number of refinement steps
NREFS = 5;

cnd = zeros(NREFS+1,length(coeffvals));

for k=1:length(coeffvals)
  fprintf('>>>> tau = %f <<<<<\n',coeffvals(k));
  extev = auxpcexpcoeff(TriangleSubDom,coeffvals(k),1,NREFS);
  cnd(:,k) = extev(:,2)./extev(:,1);
end

save 'TriangleSubDomCoeffAlphaCnd' cnd coeffvals extev NREFS 

clear cnd extev;

cnd = zeros(NREFS+1,length(coeffvals));

for k=1:length(coeffvals)
  fprintf('>>>> tau = %f <<<<<\n',coeffvals(k));
  extev = auxpcexpcoeff(TriangleSubDom,1,coeffvals(k),NREFS);
  cnd(:,k) = extev(:,2)./extev(:,1);
end

save 'TriangleSubDomCoeffBetaCnd' cnd coeffvals extev NREFS 

