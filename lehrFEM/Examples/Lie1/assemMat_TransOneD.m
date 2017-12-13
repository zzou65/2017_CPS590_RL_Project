function varargout = assemMat_TransOneD(Mesh, b_Mesh)
% assemMat_TransOneD assemble transfermatrix of one forms from finer b_Mesh to mesh
%
%   A = assemMat_TransOneD(Mesh,b_Mesh) 
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_TransOneD(MESH) .... assembles the global matrix 
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
  I=zeros(nEdges*b_Mesh.Max_Nodes,1);
  J=zeros(nEdges*b_Mesh.Max_Nodes,1);
  T=zeros(nEdges*b_Mesh.Max_Nodes,1);
  
  bd=Get_BdEdges(Mesh);
  
  loc=1;
  last=0;
  for (i = 1:nEdges)
    if (Mesh.BdFlags(i)>=0) % only for nonboundary edges
     vid=Mesh.Edges(i,:);
    
     % node 1
     if (Mesh.VBdFlags(vid(1))>=0)  % non boundary nodes
      nodes=b_Mesh.AdjNodes(vid(1),:);
      nodes=nodes(nodes~=0);
      s=size(nodes,2);
      while(nodes(1)~=nCoordinates+i)
         nodes=[nodes(2:s) nodes(1)];
      end
      eid_d=zeros(s+1,1);
      eid_d(1)=b_Mesh.Vert2Edge(nodes(1),nodes(2));
      eid_d(s+1)=b_Mesh.Vert2Edge(nodes(s),nodes(1));
      P=zeros(1,s+1);
      if(b_Mesh.Edges(eid_d(1),1)==nodes(1)),  p = 1;  else    p = -1;  end 
      P(1)=p;
      for j=2:s
         eid_d(j)=b_Mesh.Vert2Edge(nodes(j),vid(1));
         if(b_Mesh.Edges(eid_d(j),1)==i),  p = 1;  else    p = -1;  end 
         P(j)= p;
      end
      if(b_Mesh.Edges(eid_d(s+1),1)==nodes(1)),  p = 1;  else    p = -1;  end 
      P(s+1)=p;
      T(loc:(loc+s))=[-1/2 ((s/2-1):-1:(-s/2+1))/s 1/2].*P;
      %T(loc:(loc+s))=[-1/2 zeros(1, s-1) 1/2].*P;
      I(loc:(loc+s))=i*ones(1,s+1);
      J(loc:(loc+s))=eid_d;
      last=loc;
      loc=loc+s+1;
      
     % else ? boundary
     end  
     
      % node 2
     if (Mesh.VBdFlags(vid(2))>=0) % non boundary nodes
      nodes=b_Mesh.AdjNodes(vid(2),:);
      nodes=nodes(nodes~=0);
      s=size(nodes,2);
      while(nodes(1)~=nCoordinates+i)
         nodes=[nodes(2:s) nodes(1)];
      end
      eid_d=zeros(s+1,1);
      eid_d(1)=b_Mesh.Vert2Edge(nodes(1),nodes(2));
      eid_d(s+1)=b_Mesh.Vert2Edge(nodes(s),nodes(1));
      P=zeros(1,s+1);
      if(b_Mesh.Edges(eid_d(1),1)==nodes(1)),  p = 1;  else    p = -1;  end 
      P(1)=p;
      for j=2:s
         eid_d(j)=b_Mesh.Vert2Edge(nodes(j),vid(2));
         if(b_Mesh.Edges(eid_d(j),1)==i),  p = 1;  else    p = -1;  end 
         P(j)=p;
      end
      if(b_Mesh.Edges(eid_d(s+1),1)==nodes(1)),  p = 1;  else    p = -1;  end 
      P(s+1)=p;
      T(loc:(loc+s))=[-1/2 ((s/2-1):-1:(-s/2+1))/s 1/2].*P;
      %T(loc:(loc+s))=[-1/2 zeros(1, s-1) 1/2].*P;
      I(loc:(loc+s))=i*ones(1,s+1);
      J(loc:(loc+s))=eid_d;
      last=loc;
      loc=loc+s+1;
     % else ? boundary
     end
    end
     
  end
  
  % Assign output arguments
  I=I(1:last);
  J=J(1:last);
  T=T(1:last);
  if (nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = T;
  else
    varargout{1} = sparse(I,J,T,nEdges,size(b_Mesh.Edges,1));      
  end  
  
return