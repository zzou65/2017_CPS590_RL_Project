function varargout = assemMat_Mass0fD(Mesh, varargin)
% assemMat_Mass0fD assemble diagonal MASS matrix of zero forms, e.g area of
% cells
%
%   A = assemMat_0fD(Mesh,varargin) 
%   A = ASSEMMat_0fD(MESH) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_Mass0fD .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemMat_MassZeroD(Mesh);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments
  nF=size(Mesh.Elements,1);
  nC=size(Mesh.Coordinates,1);
  MD=zeros(nC,1);
  I=1:nC;
  
  for (f = 1:nF)
       
     % Vertices
     vid = Mesh.Elements(f,:);
     a1 = Mesh.Coordinates(vid(1),:);
     a2 = Mesh.Coordinates(vid(2),:);
     a3 = Mesh.Coordinates(vid(3),:);
     
     % Compute element mapping

     bK = a1;
     BK = [a2-bK; ...
        a3-bK];
     det_BK = abs(det(BK));
     inv_BK = inv(BK);
     
     MD(vid)=MD(vid)+det_BK/6*[1 1 1]';
     
  end
  
  if(nargout > 1)
    varargout{1} = MD;
    varargout{2} = I;
    varargout{3} = I;
  else
    varargout{1} = sparse(I,I,MD);      
  end  
  
  return
   
  