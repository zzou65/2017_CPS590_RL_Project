function varargout = assemMat_TransOneP(Mesh, b_Mesh)
% assemMat_TransOneP assemble transfermatrix of one forms from finer b_Mesh to mesh
%
%   A = assemMat_TransOneP(Mesh,b_Mesh) 
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_TransOneP(MESH) .... assembles the global matrix 
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

  % Initialize constants

  nCoordinates=size(Mesh.Coordinates,1);
  nEdges=size(Mesh.Edges,1);
  nElements=size(Mesh.Elements,1);
  I=zeros(18*nElements,1);
  J=zeros(18*nElements,1);
  T=zeros(18*nElements,1);
  bd=Get_BdEdges(Mesh);
  non_bd=setdiff(1:nEdges,bd);
  % inner edges will be visted twice
  InEd=ones(nEdges,1);
  InEd(non_bd)=1/2;
  loc=1:6;
  pp=zeros(1,3);
  pd=zeros(1,18);
  for (i = 1:nElements)
    % vertieces  
    vid=Mesh.Elements(i,:);
    % edges, coarse grid
    eid = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];
    % midpoints
    mid=nCoordinates+eid;
    % barycener
    bid=nCoordinates+nEdges+i;
    
    % edges, fine grid
    eid_d(1)=b_Mesh.Vert2Edge(vid(1),mid(3));
    eid_d(2)=b_Mesh.Vert2Edge(mid(3),vid(2));
    eid_d(3)=b_Mesh.Vert2Edge(vid(2),mid(1));
    eid_d(4)=b_Mesh.Vert2Edge(mid(1),vid(3));
    eid_d(5)=b_Mesh.Vert2Edge(vid(3),mid(2));
    eid_d(6)=b_Mesh.Vert2Edge(mid(2),vid(1));
    eid_d(7)=b_Mesh.Vert2Edge(vid(1),bid);
    eid_d(8)=b_Mesh.Vert2Edge(mid(3),bid);
    eid_d(9)=b_Mesh.Vert2Edge(vid(2),bid);
    eid_d(10)=b_Mesh.Vert2Edge(mid(1),bid);
    eid_d(11)=b_Mesh.Vert2Edge(vid(3),bid);
    eid_d(12)=b_Mesh.Vert2Edge(mid(2),bid);
    
    % Determine the orientation
    % coarse
    if(Mesh.Edges(eid(1),1)==vid(2)),  pp(1) = 1;  else    pp(1) = -1;  end
    if(Mesh.Edges(eid(2),1)==vid(3)),  pp(2) = 1;  else    pp(2) = -1;  end
    if(Mesh.Edges(eid(3),1)==vid(1)),  pp(3) = 1;  else    pp(3) = -1;  end
    %fine
    if(b_Mesh.Edges(eid_d(1),1)== vid(1)),  pd(1) = 1;  else    pd(1) = -1;  end
    if(b_Mesh.Edges(eid_d(2),1)== mid(3)),  pd(2) = 1;  else    pd(2) = -1;  end
    if(b_Mesh.Edges(eid_d(3),1)== vid(2)),  pd(3) = 1;  else    pd(3) = -1;  end
    if(b_Mesh.Edges(eid_d(4),1)== mid(1)),  pd(4) = 1;  else    pd(4) = -1;  end
    if(b_Mesh.Edges(eid_d(5),1)== vid(3)),  pd(5) = 1;  else    pd(5) = -1;  end
    if(b_Mesh.Edges(eid_d(6),1)== mid(2)),  pd(6) = 1;  else    pd(6) = -1;  end
    if(b_Mesh.Edges(eid_d(7),1)== vid(1)),  pd(7) = 1;  else    pd(7) = -1;  end
    if(b_Mesh.Edges(eid_d(8),1)== mid(3)),  pd(8) = 1;  else    pd(8) = -1;  end
    if(b_Mesh.Edges(eid_d(9),1)== vid(2)),  pd(9) = 1;  else    pd(9) = -1;  end
    if(b_Mesh.Edges(eid_d(10),1)== mid(1)),  pd(10) = 1;  else    pd(10) = -1;  end
    if(b_Mesh.Edges(eid_d(11),1)== vid(3)),  pd(11) = 1;  else    pd(11) = -1;  end
    if(b_Mesh.Edges(eid_d(12),1)== mid(2)),  pd(12) = 1;  else    pd(12) = -1;  end
    % Edge 1
    I(loc)=eid(1)*ones(1,6);
    J(loc)=[eid_d(3) eid_d(4) eid_d(9) eid_d(11) eid_d(8) eid_d(12)];
    T(loc)=pp(1)*[InEd(eid(1))*1/2*pd(3) InEd(eid(1))*1/2*pd(4) 1/3*pd(9) -1/3*pd(11) 1/6*pd(8) -1/6*pd(12)];
    loc=loc+6;
    % Edge 2
    I(loc)=eid(2)*ones(1,6);
    J(loc)=[eid_d(5) eid_d(6) eid_d(11) eid_d(7) eid_d(10) eid_d(8)];
    T(loc)=pp(2)*[InEd(eid(2))*1/2*pd(5) InEd(eid(2))*1/2*pd(6) 1/3*pd(11) -1/3*pd(7) 1/6*pd(10) -1/6*pd(8)];
    loc=loc+6;
    % Edge 3
    I(loc)=eid(3)*ones(1,6);
    J(loc)=[eid_d(1) eid_d(2) eid_d(7) eid_d(9) eid_d(12) eid_d(10)];
    T(loc)=pp(3)*[InEd(eid(3))*1/2*pd(1) InEd(eid(3))*1/2*pd(2) 1/3*pd(7) -1/3*pd(9) 1/6*pd(12) -1/6*pd(10)];
    loc=loc+6;
  end
  
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = T;
  else
    varargout{1} = sparse(I,J,T);      
  end  
  
return 