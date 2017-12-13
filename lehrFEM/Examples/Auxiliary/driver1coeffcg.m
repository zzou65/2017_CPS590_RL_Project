% Numerical experiment I: driver routine

clear all;

% Scaling parameter
coeffvals = [0.001,0.01,0.1,2,5,10,20,50,100,200,1000];

% Load mesh
load 'TriangleSubDom.mat';

% Number of refinement steps
NREFS = 5;

itnum = zeros(NREFS+1,length(coeffvals));

for k=1:length(coeffvals)
  TriangleSubDom.cg = true;
  fprintf('>>>> tau = %f <<<<<\n',coeffvals(k));
  extev = auxpcexpcoeff(TriangleSubDom,coeffvals(k),1,NREFS);
  itnum(:,k) = extev(:,1);
end

save 'TriangleSubDomCoeffAlphaCG' itnum coeffvals extev NREFS 

clear itnum extev;

itnum = zeros(NREFS+1,length(coeffvals));

for k=1:length(coeffvals)
  TriangleSubDom.cg = true;
  fprintf('>>>> tau = %f <<<<<\n',coeffvals(k));
  extev = auxpcexpcoeff(TriangleSubDom,1,coeffvals(k),NREFS);
  itnum(:,k) = extev(:,1);
end

save 'TriangleSubDomCoeffBetaCG' itnum coeffvals extev NREFS 
