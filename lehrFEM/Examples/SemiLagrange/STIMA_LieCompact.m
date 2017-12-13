function Aloc = STIMA_LieCompact(Vertices,ElemInfo,V_HANDLE,QuadRule,varargin)
% STIMA_LieCompact computes element contribution of modified Lie-derivative 
%
%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

Aloc = zeros(3,3);
% Compute element mapping
P1 = Vertices(1,:);
P2 = Vertices(2,:);
P3 = Vertices(3,:);

BK = [ P2 - P1 ; P3 - P1 ]; % transpose of transformation matrix
det_BK = abs(det(BK));     % twice the area of the triagle

% Evaluate coefficient function at vertices
Fval = V_HANDLE(Vertices,ElemInfo,varargin{:});

% calculate barycentric coordinates of velocity
Fval(1,:)= (P1+Fval(1,:)-P1) * inv(BK);
Fval(2,:)= (P2+Fval(2,:)-P1) * inv(BK);
Fval(3,:)= (P3+Fval(3,:)-P1) * inv(BK);

Fval=[1-sum(Fval,2),Fval];

% edge [2 3]
if (Fval(2,1)<0)   Aloc(1,:)=[Fval(2,2)-1 -Fval(2,1) 0]; end;
if (Fval(2,1)==0)   Aloc(1,:)=1/2*[Fval(2,2)-1 -Fval(2,1) 0]; end;
if (Fval(3,1)<0)   Aloc(1,:)=Aloc(1,:)+[Fval(3,3)-1 0 -Fval(3,1)]; end;
if (Fval(3,1)==0)   Aloc(1,:)=Aloc(1,:)+1/2*[Fval(3,3)-1 0 -Fval(3,1)]; end;
% edge [3 1]
if (Fval(3,2)<0)   Aloc(2,:)=[0 Fval(3,3)-1 -Fval(3,2)]; end;
if (Fval(3,2)==0)   Aloc(2,:)=1/2*[0 Fval(3,3)-1 -Fval(3,2)]; end;
if (Fval(1,2)<0)   Aloc(2,:)=Aloc(2,:)+[-Fval(1,2) Fval(1,1)-1 0]; end;
if (Fval(1,2)==0)   Aloc(2,:)=Aloc(2,:)+1/2*[-Fval(1,2) Fval(1,1)-1 0]; end;
% edge [1 2]
if (Fval(1,3)<0)   Aloc(3,:)=[-Fval(1,3) 0 Fval(1,1)-1]; end;
if (Fval(1,3)==0)   Aloc(3,:)=1/2*[-Fval(1,3) 0 Fval(1,1)-1]; end;
if (Fval(2,3)<0)   Aloc(3,:)=Aloc(3,:)+[0 -Fval(2,3) Fval(2,2)-1]; end;
if (Fval(2,3)==0)   Aloc(3,:)=Aloc(3,:)+1/2*[0 -Fval(2,3) Fval(2,2)-1]; end;

return