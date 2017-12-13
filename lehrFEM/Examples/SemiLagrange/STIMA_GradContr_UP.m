function B = STIMA_GradContr_UP(Vertices,ElemInfo,V_HANDLE,QuadRule,varargin)
%% bloedsinn

% STIMA_GradContr_UP computes element contribution of grad v u term in
% the sense of cochains, e.g.  mapping of edge DOFS to edge DOFS in
% evaluating  the intergralvalues \int_{edge_1} grad v b_edge_j terms.
% Here we use the trace from the upwind direction.
%
% !!!!! This can not be assembled element by element !!!!! Hence this
% implementation is wrong !!!!
%
%   MLOC = STIMA_ContrRot_UP(VERTICES ...) computes element contribution of v x curl u term matrix using 
%   Whitney 1-forms finite elements.

%   Copyright 2009-2009 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

B=zeros(3,3);
H=zeros(3,3);

P1 = Vertices(1,:);
P2 = Vertices(2,:);
P3 = Vertices(3,:);
% Compute element mapping

bK = P1;
BK = [P2-bK; ...
  P3-bK];
det_BK = abs(det(BK));
inv_BK = inv(BK);

% Extract nodal vectors
v1=-V_HANDLE(P1);
v2=-V_HANDLE(P2);
v3=-V_HANDLE(P3);

% Get barycentric coordinates of quadrature points
baryc= [1 0 0;  0 1 0; 0 0 1];

% Compute constant gradients of barycentric coordinate functions
g1 = [P2(2)-P3(2);P3(1)-P2(1)]/det_BK;
g2 = [P3(2)-P1(2);P1(1)-P3(1)]/det_BK;
g3 = [P1(2)-P2(2);P2(1)-P1(1)]/det_BK;

% Evaluate basis functions at quadrature points
% the rows of b(i,:) store the value of the the i-th
% basis function at the quadrature points
b1 = baryc(:,2)*g3'-baryc(:,3)*g2';
b2 = baryc(:,3)*g1'-baryc(:,1)*g3';
b3 = baryc(:,1)*g2'-baryc(:,2)*g1';

% Compute decision variables Theta(T,

% for first vertex
y = P1 + v1;
yhat = (y-bK)*inv_BK;
if(yhat(1) >= 0 && yhat(2) >= 0)
    H(1,:)=[0 0 0];
    H(2,:)=[v1*b1(1,:)' v1*b2(1,:)' v1*b3(1,:)'];
    H(3,:)=-H(2,:);
    if (yhat(1)<eps || yhat(2)<eps)
        H=1/2 * H;
    end
    B=B+H;
end

% for second vertex
y = P2 + v2;
yhat = (y-bK)*inv_BK;
if(yhat(2) >= 0 && sum(yhat) <= 1)
    H(1,:)=-[v2*b1(2,:)' v2*b2(2,:)' v2*b3(2,:)'];
    H(2,:)=[0 0 0]; 
    H(3,:)=-H(1,:);
       if (yhat(2)<eps || abs(sum(yhat)-1) <eps )
           H=1/2 * H;
    end
    B=B+H;
end

% for third index
y = P3 + v3;
yhat = (y-bK)*inv_BK;
if(yhat(1) >= 0 && sum(yhat) <= 1)
    H(1,:)=[v3*b1(3,:)' v3*b2(3,:)' v3*b3(3,:)'];
    H(2,:)=-H(1,:);
    H(3,:)=[0 0 0];
    if (yhat(1)<eps || abs(sum(yhat)-1) < eps)
        H=1/2 * H;
    end
    B=B+H;
end
B=-B;