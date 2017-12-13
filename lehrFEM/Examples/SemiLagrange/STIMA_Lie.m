function Mloc = STIMA_Lie(Vertices,ElemInfo,V_HANDLE,JAC_HANDLE,varargin)
% STIMA_Lie computes element contribution of -v x curl u +grad v u  term in
% using the identity -v x rot u + grad v u = Dv u +D u v + u x curl v 
%
%   Copyright 2009 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping

  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ]; % transpose of transformation matrix
  det_BK = abs(det(BK));     % twice the area of the triagle
  
  % Compute constant gradients of barycentric coordinate functions
  g1 = [P2(2)-P3(2);P3(1)-P2(1)]/det_BK;
  g2 = [P3(2)-P1(2);P1(1)-P3(1)]/det_BK;
  g3 = [P1(2)-P2(2);P2(1)-P1(1)]/det_BK;
  
  % Get barycentric coordinates of quadrature points 
  baryc= [1/3 1/3 1/3];

  % Quadrature points in actual element
  % stored as rows of a matrix
  bary = 1/3*(P1+P2+P3);

  % Evaluate coefficient function at quadrature nodes
  Velo = V_HANDLE(bary,ElemInfo,varargin{:});
  JVelo = JAC_HANDLE(bary,ElemInfo,varargin{:});
  JVelo=reshape(JVelo,2,2);

  % Evaluate basis functions at quadrature points
  % the rows of b(i) store the value of the the i-th
  % basis function at the quadrature points
  b1 = baryc(:,2)*g3'-baryc(:,3)*g2';
  b2 = baryc(:,3)*g1'-baryc(:,1)*g3';
  b3 = baryc(:,1)*g2'-baryc(:,2)*g1';
    
  % Compute local matrix
  % (Dv u, u')
  Mloc1(:,:) = det_BK/2*...
      [b1 *JVelo*b1' b2*JVelo*b1' b3 *JVelo*b1';...
       b1 *JVelo*b2' b2*JVelo*b2' b3 *JVelo*b2';...
       b1 *JVelo*b3' b2*JVelo*b3' b3 *JVelo*b3']';
   
  % (Du v, u')=-0.5*(v x curl u,u')
  Ju=[0 1; -1 0];
  Mloc2(:,:) = 1/2*...
      [Velo*Ju*b1' Velo*Ju*b1' Velo*Ju*b1';...
       Velo*Ju*b2' Velo*Ju*b2' Velo*Ju*b2';...
       Velo*Ju*b3' Velo*Ju*b3' Velo*Ju*b3'];

 % (u x curl v, u') = (u' x u, curl v)
 curlv=JVelo(2,1)-JVelo(1,2);
 Ju=[ 0 -1; 1 0];
 Mloc3(:,:) = curlv*det_BK/2*...
     [b1*Ju*b1' b2*Ju*b1' b3*Ju*b1';...
     b1*Ju*b2' b2*Ju*b2' b3*Ju*b2';...
     b1*Ju*b3' b2*Ju*b3' b3*Ju*b3'];
 
 Mloc = Mloc1+Mloc2+Mloc3;

 return