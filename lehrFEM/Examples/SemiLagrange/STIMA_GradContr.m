function Mloc = STIMA_GradContr(Vertices,ElemInfo,V_HANDLE,QuadRule,varargin)
%% bloedsinn

% STIMA_GradContr computes element entries for (grad v u, u') in DG-like
% setting, (v u, div u')_T vanishes and only boundary terms remain.
%
%   MLOC = STIMA_GradContr(VERTICES,...) computes ... matrix using 
%   Whitney 1-forms finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
% 
%   ElemInfo (not used)
%
%   V_HANDLE 
%   Example:
%

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%   
  Mloc=zeros(3,3); 

  % Compute element mapping
   
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ];   % transpose of transformation matrix
  det_BK = abs(det(BK));        % twice the area of the triagle  
  inv_BK = inv(BK); 
  
  % Compute constant gradients of barycentric coordinate functions

  g1 = [P2(2)-P3(2);P3(1)-P2(1)]/det_BK;
  g2 = [P3(2)-P1(2);P1(1)-P3(1)]/det_BK;
  g3 = [P1(2)-P2(2);P2(1)-P1(1)]/det_BK;
    
  % Quadrature points in actual element
  % stored as rows of a matrix
  
  nPoints = size(QuadRule.w,1);
  
  x3 = QuadRule.x*(P2-P1) + ones(nPoints,1)*P1;
  x1 = QuadRule.x*(P3-P2) + ones(nPoints,1)*P2;
  x2 = QuadRule.x*(P1-P3) + ones(nPoints,1)*P3;
   
  % Get barycentric coordinates of quadrature points 
  
  %baryc= [[x1 ; x2 ; x3],1-sum([x1 ; x2 ; x3],2)];
  baryc= [zeros(nPoints,1) 1-QuadRule.x QuadRule.x;...
              QuadRule.x zeros(nPoints,1) 1-QuadRule.x;...
              1-QuadRule.x QuadRule.x zeros(nPoints,1)];

  % Evaluate coefficient function at quadrature nodes
  Fval3 = V_HANDLE(x3,ElemInfo,varargin{:});
  Fval1 = V_HANDLE(x1,ElemInfo,varargin{:});
  Fval2 = V_HANDLE(x2,ElemInfo,varargin{:});
  Fval=[Fval1; Fval2; Fval3];
  
  % Evaluate basis functions at quadrature points
  % the rows of b(i) store the value of the the i-th
  % basis function at the quadrature points
  
  b1 = baryc(:,2)*g3'-baryc(:,3)*g2';
  b2 = baryc(:,3)*g1'-baryc(:,1)*g3';
  b3 = baryc(:,1)*g2'-baryc(:,2)*g1';

  % edge length
  dS=[norm(P3-P2),norm(P1-P3),norm(P2-P1)];
  
  % Quadrature weights
  w=[dS(1)*QuadRule.w ; dS(2)*QuadRule.w; dS(3)*QuadRule.w]; 
  
  % Compute outward normals
  N1=[1,1]*inv_BK';
  N1=N1/norm(N1);
  N2=[-1,0]*inv_BK';
  N2=N2/norm(N2);
  N3=[0,-1]*inv_BK';
  N3=N3/norm(N3);

  % compute auxiliary products
  %   
  %   b1n1=sum(b1(1:nPoints).*(ones(nPoints,1)*N1),2);
  %   b1n2=sum(b1(nPoints+1:nPoints).*(ones(nPoints,1)*N2),2); 
  %   b1n3=sum(b1(2*nPoints+1:nPoints).*(ones(nPoints,1)*N3),2);
  %   
  %   b2n1=sum(b2(1:nPoints).*(ones(nPoints,1)*N1),2);
  %   b2n2=sum(b2(nPoints+1:nPoints).*(ones(nPoints,1)*N2),2);
  %   b2n3=sum(b2(2*nPoints+1:nPoints).*(ones(nPoints,1)*N3),2);
  %   
  %   b3n1=sum(b3(1:nPoints).*(ones(nPoints,1)*N1),2);
  %   b3n2=sum(b3(nPoints+1:nPoints).*(ones(nPoints,1)*N2),2);
  %   b3n3=sum(b3(2*nPoints+1:nPoints).*(ones(nPoints,1)*N3),2);
  %   
  b1n=sum(b1.*([ones(nPoints,1)*N1;...
                        ones(nPoints,1)*N2;...
                        ones(nPoints,1)*N3]),2);
 
  b2n=sum(b2.*([ones(nPoints,1)*N1;...
                        ones(nPoints,1)*N2;...
                        ones(nPoints,1)*N3]),2);
 
  b3n=sum(b3.*([ones(nPoints,1)*N1;...
                        ones(nPoints,1)*N2;...
                        ones(nPoints,1)*N3]),2);
  
  vb1=sum(Fval.*b1,2);
  vb2=sum(Fval.*b2,2);
  vb3=sum(Fval.*b3,2);
  
  Mloc(1,1)=sum(b1n.*vb1.*w);
  Mloc(2,1)=sum(b2n.*vb1.*w);
  Mloc(3,1)=sum(b3n.*vb1.*w);
  
  Mloc(1,2)=sum(b1n.*vb2.*w);
  Mloc(2,2)=sum(b2n.*vb2.*w);
  Mloc(3,2)=sum(b3n.*vb2.*w);
  
  Mloc(1,3)=sum(b1n.*vb3.*w);
  Mloc(2,3)=sum(b2n.*vb3.*w);
  Mloc(3,3)=sum(b3n.*vb3.*w);
  Mloc=Mloc;
  return;
  