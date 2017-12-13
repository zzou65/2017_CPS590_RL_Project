function L = assemCochain_2f(Mesh, FHandle, QuadRule)
% assemCochain_1f assemble Cochain of 2 forms
%
%   L = ASSEMCochain_1f(MESH,FHandel) .... 
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemMat_MassTwoD(Mesh);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments
 
 nElements=size(Mesh.Elements,1);
 L = zeros(nElements,1);
 nPts = size(QuadRule.w,1);
 
 % Assemble element contributions
  
 for i = 1:nElements

        vid = Mesh.Elements(i,:);
        vert=Mesh.Coordinates(vid,:);
        
        %  element mapping 
        
        bK=vert(1,:);
        BK=[vert(2,:)-vert(1,:);vert(3,:)-vert(1,:)];
        det_BK=abs(det(BK));
        
        %Quadrature Points
        
        x = QuadRule.x*BK + ones(nPts,1)*bK;
    
        % Compute load data
    
        FVal = FHandle(x);
 
       % Add contributions to global load vector
        
        L(i) =  sum(QuadRule.w.*FVal)*det_BK;
end
  
return
   
  