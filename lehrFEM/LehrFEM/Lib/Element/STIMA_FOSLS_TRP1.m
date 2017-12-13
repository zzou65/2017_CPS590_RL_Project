function Aloc = STIMA_FOSLS_TRP1(Vertices,ElemInfo,QuadRule,varargin)
% STIMA_FOSLS_TRP1 Element stiffness matrix for the laplace problem in FOSLS formulation.
%
%   ALOC = STIMA_FOSLS_TRP1(VERTICES,ELEMINFO,QUADRULE) computes the
%   element stiffness matrix for the laplace problem in FOSLS formulation. 
%
%   VERTICES is a 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%   
%   Example:
%
%   Aloc = STIMA_FOSLS_TRP1([0 0; 1 0; 0 1],0,P7O6());
%

%   Copyright 2005-2006 Patrick Meury & Kah Ling Sia & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  MU_Handle = @(x,varargin)ones(size(x,1),1);
  Aloc = zeros(6,6);
  
  % Compute element mapping
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  det_BK = abs(det(BK));
   
  % Compute constant gradients of barycentric coordinate functions
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  
  a1=P1;
  a2=P2;
  a3=P3;

  
 P=[ 2*a1(2)-a3(2)-a2(2)  a2(2)-a3(2)          a3(2)-a2(2);...
   -2*a1(1)+a3(1)+a2(1) -a2(1)+a3(1)         -a3(1)+a2(1); ...
    a1(2)-a3(2)          2*a2(2)-a3(2)-a1(2)  a3(2)-a1(2); ...
   -a1(1)+a3(1)         -2*a2(1)+a3(1)+a1(1) -a3(1)+a1(1); ...
    a1(2)-a2(2)          a2(2)-a1(2)          2*a3(2)-a1(2)-a2(2);... 
   -a1(1)+a2(1)         -a2(1)+a1(1)         -2*a3(1)+a1(1)+a2(1);
   ];


  
  g1 = [P2(2)-P3(2);P3(1)-P2(1)]/det_BK;
  g2 = [P3(2)-P1(2);P1(1)-P3(1)]/det_BK;
  g3 = [P1(2)-P2(2);P2(1)-P1(1)]/det_BK;
  
  %barycenter
  
  baryc=1/3*(P1+P2+P3);
  
  % Compute element stiffness matrix
  
  % mass term of flux discretisation
  Aloc(1:3,1:3)=MASS_W1F(Vertices,0,MU_Handle,QuadRule,varargin);
%  Aloc(1:3,1:3)=1/(24*det_BK)*P'*P;
  % (div , div )-term
  Aloc(1:3,1:3)=Aloc(1:3,1:3)+2/(det_BK)*ones(3,3);
  
  % laplaceian term
  Aloc(4,4)=det_BK/2.*g1'*g1;
  Aloc(4,5)=det_BK/2.*g1'*g2;
  Aloc(4,6)=det_BK/2.*g1'*g3;
  Aloc(5,4)=det_BK/2.*g2'*g1;
  Aloc(5,5)=det_BK/2.*g2'*g2;
  Aloc(5,6)=det_BK/2.*g2'*g3;
  Aloc(6,4)=det_BK/2.*g3'*g1;
  Aloc(6,5)=det_BK/2.*g3'*g2;
  Aloc(6,6)=det_BK/2.*g3'*g3;
  
  % upper off diagonal block
  Aloc(1,4)=(baryc-P1)*g1/2;
  Aloc(1,5)=(baryc-P1)*g2/2;
  Aloc(1,6)=(baryc-P1)*g3/2;
  Aloc(2,4)=(baryc-P2)*g1/2;
  Aloc(2,5)=(baryc-P2)*g2/2;
  Aloc(2,6)=(baryc-P2)*g3/2;
  Aloc(3,4)=(baryc-P3)*g1/2;
  Aloc(3,5)=(baryc-P3)*g2/2;
  Aloc(3,6)=(baryc-P3)*g3/2;
  
  % lower off diagonal block
  Aloc(4:6,1:3)=Aloc(1:3,4:6)';
  Aloc=Aloc;
  
return