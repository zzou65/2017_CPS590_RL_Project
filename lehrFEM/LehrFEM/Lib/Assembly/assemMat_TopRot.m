function varargout = assemMat_TopRot(Mesh, varargin)
% assemMat_TopGrad Assemble topological rotation.
%
%   A = assemMat_TopGrad(Mesh,varargin) 
%   A = ASSEMMat_TopGrad(MESH) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_TopGrad(MESH) .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemMat_TopRot(Mesh);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments

  nF=size(Mesh.Elements,1);
  I1=zeros(nF,1);
  I2=zeros(nF,1);
  I3=zeros(nF,1);
  S1=zeros(nF,1);
  S2=zeros(nF,1);
  S3=zeros(nF,1);
  
  for (f = 1:nF)
     
    vidx = Mesh.Elements(f,:);
    Vertices = Mesh.Coordinates(vidx,:);
      
    % Extract global edge numbers
    
    I1(f) = Mesh.Vert2Edge(Mesh.Elements(f,2),Mesh.Elements(f,3));
    I2(f) = Mesh.Vert2Edge(Mesh.Elements(f,3),Mesh.Elements(f,1));
    I3(f) = Mesh.Vert2Edge(Mesh.Elements(f,1),Mesh.Elements(f,2));
       
    % Determine the orientation
    
    if (Mesh.Edges(I1(f),1)==vidx(2)),  S1(f) = 1;  else   S1(f) = -1;  end
    if (Mesh.Edges(I2(f),1)==vidx(3)),  S2(f) = 1;  else   S2(f) = -1;  end
    if (Mesh.Edges(I3(f),1)==vidx(1)),  S3(f) = 1;  else   S3(f) = -1;  end
    
  end
  
  if(nargout > 1)
    varargout{1} = [1:nF; 1:nF; 1:nF]';
    varargout{2} = [I1;I2;I3];
    varargout{3} = [S1;S2:S3];
  else
    varargout{1} = sparse([1:nF 1:nF 1:nF]',[I1;I2;I3],[S1;S2;S3]);      
  end  
  
  return
   
  
