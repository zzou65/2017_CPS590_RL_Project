function [ContrRot C]=assemData_ContrRot_UpQuad(Mesh, vHandle, varargin)
%  assembles data for consistent scheme for -v x curl u based on upwind quadrature using P3O3() 
%  see assemMat_ContrRot_UPQuad(

%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  TOL=eps;
  nCoordinates = size(Mesh.Coordinates,1);
  nEdges =size(Mesh.Edges,1);
  nElements=size(Mesh.Elements,1);
 
  % Preallocate memory
  %I = zeros(9*nElements,1);
  %J = zeros(9*nElements,1);
  A = zeros(nEdges,3);
  C = zeros(nEdges,3);
  
  % look for upwind nodes on each triangle
  rot=[0 1; -1 0];
  loc=1:9;
  
  for i = 1:nElements
      
     B=zeros(3,3);
     
     % Vertices
     vid = Mesh.Elements(i,:);
     a1 = Mesh.Coordinates(vid(1),:);
     a2 = Mesh.Coordinates(vid(2),:);
     a3 = Mesh.Coordinates(vid(3),:);
     
     % rotated egdes = inward normal
     n1=(a3-a2)*rot;
     n2=(a1-a3)*rot;
     n3=(a2-a1)*rot;
     
     % Extract global edge numbers
     eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];
        
    if (Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else   p1 = -1;  end
    if (Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else   p2 = -1;  end
    if (Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else   p3 = -1;  end       
       
     % Compute element mapping

     bK = a1;
     BK = [a2-bK; ...
        a3-bK];
     det_BK = abs(det(BK));
  
     % Extract nodal vectors
     v1 =-vHandle((a2+a3)/2);
     v2 =-vHandle((a3+a1)/2);
     v3 =-vHandle((a1+a2)/2);

     % Compute decision variables Theta
     
     % for first edge
     if ((v1)*n1'>TOL*norm(n1)*norm(v1))
         A(eidx(1),:)=[p1*2/det_BK p2*2/det_BK p3*2/det_BK];
         C(eidx(1),:)=eidx;
     end
     
     % for second edge
     if ((v2)*n2'>TOL*norm(v2)*norm(n2))
         A(eidx(2),:)=[p1*2/det_BK p2*2/det_BK p3*2/det_BK];
         C(eidx(2),:)=eidx;
     end
     
     % for third edge
     if ((v3)*n3'>TOL*norm(v3)*norm(n3))
         A(eidx(3),:)=[p1*2/det_BK p2*2/det_BK p3*2/det_BK];
         C(eidx(3),:)=eidx;
     end
     
  end
  ContrRot=A;
%   end
  
return
