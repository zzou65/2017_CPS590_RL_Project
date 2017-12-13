function wp= assemWedge_1f1f(Mesh,E1,E2)
% assem assemWedge_1f1f assemble transfermatrix of one forms from finer b_Mesh to mesh
%
%   A = assemWedge_1f1f(Mesh,E1,E2) 
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEM assemWedge_1f1f(MESH) .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   b_Mesh =refin_BAR(Mesh);
%   A = assemMat_TransOneP(Mesh,b_Mesh);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

nCoordinates=size(Mesh.Coordinates,1);
wp=zeros(nCoordinates,1);
for i= 1: nCoordinates
   AdjNodes=Mesh.AdjNodes(i,:);
   AdjNodes=AdjNodes(AdjNodes~=0);
   nAdjNodes=size(AdjNodes,2);
   weight=0;
   for j=AdjNodes
       edge=Mesh.Vert2Edge(i,j);
       Elem=Mesh.Edge2Elem(edge,:);
       if (Elem(1)~=0)
           vid=Mesh.Elements(Elem(1),:);
           k=setdiff(vid,[i,j]);
           p=Mesh.Coordinates(vid,:);
           det_Bk=abs(det([p(1,:)-p(2,:);p(1,:)-p(3,:)]));
           f=Mesh.Vert2Edge(i,k);
           
           % Handle boundary edges
           a=1;
           b=1;
           if (Mesh.BdFlags(edge)<0)
               a=2;
           end
           if (Mesh.BdFlags(f)<0)
               b=2;
           end
           % Determine the orientation
    
           if(Mesh.Edges(edge,2)==Mesh.Edges(f,1) || Mesh.Edges(edge,1)==Mesh.Edges(f,2)),  p1 = -1;  else    p1 = 1;  end
           
           % weight of primal and dual Edge element Basis functions on
           % Intersection with cell attached to vertice i
%           weight=weight+a/(6*det_Bk)*E1(edge)*E2(edge)+p1*b/(3*det_Bk)*E1(edge)*E2(f);
           
           [m1 m2]=WEIGHT_WEDGE([Mesh.Coordinates(i,:); ...
               Mesh.Coordinates(j,:); ...
               Mesh.Coordinates(k,:)]);
           weight=weight+a*m1*E1(edge)*E2(edge)+p1*b*m2*E1(edge)*E2(f); 
           %weight=weight+1/(72)*E1(edge)*E2(edge)+p1/(36)*E1(edge)*E2(f);
%            weight=weight+1/(3*det_Bk)*E1(edge)*E2(edge);
       end
       if (Elem(2)~=0)
           vid=Mesh.Elements(Elem(2),:);
           k=setdiff(vid,[i,j]);
           det_Bk=abs(det([p(1,:)-p(2,:);p(1,:)-p(3,:)]));
           f=Mesh.Vert2Edge(i,k);

           % Handle boundary edges
           a=1;
           b=1;
           if (Mesh.BdFlags(edge)<0)
               a=2;
           end
           if (Mesh.BdFlags(f)<0)
               b=2;
           end
           % Determine the orientation
    
           if(Mesh.Edges(edge,2)==Mesh.Edges(f,1) || Mesh.Edges(edge,1)==Mesh.Edges(f,2)),  p1 = -1;  else    p1 = 1;  end
           
           % weight of primal and dual Edge element Basis functions on
           % Intersection with cell attached to vertice i
           
           %weight=weight+a/(6*det_Bk)*E1(edge)*E2(edge)+p1*b/(3*det_Bk)*E1(edge)*E2(f);
           
           [m1 m2]=WEIGHT_WEDGE([Mesh.Coordinates(i,:); ...
               Mesh.Coordinates(j,:); ...
               Mesh.Coordinates(k,:)]);
           weight=weight+a*m1*E1(edge)*E2(edge)+p1*b*m2*E1(edge)*E2(f); 
           
           %weight=weight+1/(72)*E1(edge)*E2(edge)+p1/(36)*E1(edge)*E2(f);
%            weight=weight+1/(3*det_Bk)*E1(edge)*E2(edge);
       end
       wp(i)=weight;
   end
end
 