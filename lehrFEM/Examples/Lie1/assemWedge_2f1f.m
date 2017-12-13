function wp= assemWedge_2f1f(Mesh,B,E)
% assem assemWedge_2f1f assemble wedge product of two one forms
%
%   wp = assemWedge_2f1f(Mesh,B,E2) returns n-1 array wp which is the
%   discrte wedge product of p-1 array B and m-1 array E2. n number of 
%   vertices, m number of edges, p number of elements
%
%   Example:
%
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


nCoordinates=size(Mesh.Coordinates,1);
nElem=size(Mesh.Elements,1);
wp=zeros(nCoordinates,2);
for i= 1: nElem
  
   %vertices
   
   vid=Mesh.Elements(i,:);
   Vertices=Mesh.Coordinates(vid,:);
   
   % Compute element contributions corresponding to the vertices
    
   [Aloc1 Aloc2 Aloc3] = WEIGHT_WEDGE2f1f(Vertices);
   
   % Extract global edge numbers
    
   eidx = [Mesh.Vert2Edge(vid(2),vid(3)) ...
           Mesh.Vert2Edge(vid(3),vid(1)) ...
           Mesh.Vert2Edge(vid(1),vid(2))];
       
    % Determine the orientation
    
    if(Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else    p1 = -1;  end
    if(Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else    p3 = -1;  end
    
    Peori = [p1; p2 ;p3]*ones(1,2); % scaling matrix taking into account orientations
    
    Aloc1 = Peori.*Aloc1;
    Aloc2 = Peori.*Aloc2;
    Aloc3 = Peori.*Aloc3;
    
    % Add contributions to stiffness matrix
    
    wp(vid(1),:)=wp(vid(1),:)+B(i)*sum((E(eidx))*ones(1,2).*Aloc1,1);
    wp(vid(2),:)=wp(vid(2),:)+B(i)*sum((E(eidx))*ones(1,2).*Aloc2,1);
    wp(vid(3),:)=wp(vid(3),:)+B(i)*sum((E(eidx))*ones(1,2).*Aloc3,1);
  end
end
 