function varargout = assemMat_MassOneD(Mesh, varargin)
% assemMat_MassOneD Assemble diagonal MASS matrix.
%
%   A = assemMat_MassOneD(Mesh,varargin) 
%   A = ASSEMMat_MassOneD(MESH) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_MassOneD
%    D .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemMat_MassOneD(Mesh);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants

% Assign output arguments
  nE=size(Mesh.Edges,1);
  MD=zeros(nE,1);
  I=1:nE;
  
  for (e = 1:nE)
     
     %Elements
     
     elem=Mesh.Edge2Elem(e,:);
     
     % vertices of dual edges 
     % inner edges 
     if ( elem(1) ~= 0 &&  elem(2) ~=0)
     
         bary1=1/3*sum(Mesh.Coordinates(Mesh.Elements(elem(1),:),:),1);
         
         bary2=1/3*sum(Mesh.Coordinates(Mesh.Elements(elem(2),:),:),1);
     end
     % boundary edges
     if ( elem(1) == 0)
     
         bary1=1/2*sum(Mesh.Coordinates(Mesh.Edges(e,:),:),1);
         
         bary2=1/3*sum(Mesh.Coordinates(Mesh.Elements(elem(2),:),:),1);
     end

     if (elem(2) ==0)
     
         bary1=1/3*sum(Mesh.Coordinates(Mesh.Elements(elem(1),:),:),1);
         
         bary2=1/2*sum(Mesh.Coordinates(Mesh.Edges(e,:),:),1);
        
     end
     
     MD(e)=norm(bary1-bary2)/(norm(Mesh.Coordinates(Mesh.Edges(e,1),:)-Mesh.Coordinates(Mesh.Edges(e,2),:)));
  %   MD(e)=1/MD(e);
     
  end
  
  if(nargout > 1)
    varargout{1} = MD;
    varargout{2} = I;
    varargout{3} = I;
  else
    varargout{1} = sparse(I,I,MD,nE,nE);      
  end
  return
   
  