function wp= assemWedge_1f1f(Mesh,E1,E2)
% assem assemWedge_1f1f assemble wedge product of two one forms
%
%   wp = assemWedge_1f1f(Mesh,E1,E2) returns n-1 array wp which is the
%   discrte wedge product of m-1 arrays E1 and E2. n number of vertices, m
%   number od edges
%
%   Example:
%
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

nCoordinates=size(Mesh.Coordinates,1);
nElem=size(Mesh.Elements,1);
wp=zeros(nCoordinates,1);
for i= 1: nElem
  
   %vertices
   
   vid=Mesh.Elements(i,:);
   Vertices=Mesh.Coordinates(vid,:);
   
   % Compute element contributions corresponding to the edges
    
   [Aloc1 Aloc2 Aloc3] = WEIGHT_WEDGE(Vertices);
   %[Aloc1 Aloc2 Aloc3] = WEIGHT_WEDGE2(Vertices);
   % Extract global edge numbers
    
   eidx = [Mesh.Vert2Edge(vid(2),vid(3)) ...
           Mesh.Vert2Edge(vid(3),vid(1)) ...
           Mesh.Vert2Edge(vid(1),vid(2))];
       
    % Determine the orientation
    
    if(Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else    p1 = -1;  end
    if(Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else    p3 = -1;  end
    
    Peori = diag([p1 p2 p3]); % scaling matrix taking into account orientations
    Aloc1 = Peori*Aloc1*Peori;
    Aloc2 = Peori*Aloc2*Peori;
    Aloc3 = Peori*Aloc3*Peori;
    
    % Add contributions to stiffness matrix
%     
%     wp(vid(1))=wp(vid(1))+E1(eidx)'*Aloc1*E2(eidx)*sigma(vid(1));
%     wp(vid(2))=wp(vid(2))+E1(eidx)'*Aloc2*E2(eidx)*sigma(vid(2));
%     wp(vid(3))=wp(vid(3))+E1(eidx)'*Aloc3*E2(eidx)*sigma(vid(3));
    wp(vid(1))=wp(vid(1))+E1(eidx)'*Aloc1*E2(eidx);
    wp(vid(2))=wp(vid(2))+E1(eidx)'*Aloc2*E2(eidx);
    wp(vid(3))=wp(vid(3))+E1(eidx)'*Aloc3*E2(eidx);
  end
end
 