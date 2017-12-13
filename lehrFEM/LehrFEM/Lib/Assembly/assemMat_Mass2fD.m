function varargout = assemMat_Mass1fD(Mesh, varargin)
% assemMat_MassTwoD Assemble diagonal MASS matrix of 2 forms
%
%   A = assemMat_Mass2fD(Mesh,varargin) 
%   A = ASSEMMat_Mass2fD(MESH) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_Mass2f
%    D .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemMat_Mass2fD(Mesh);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments
  nF=size(Mesh.Elements,1);
  nC=size(Mesh.Coordinates,1);
  MD=zeros(nF,1);
  I=1:nF;
  
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
     
     MD(f)=2/det_BK;
     
  end
  
  if(nargout > 1)
    varargout{1} = MD;
    varargout{2} = I;
    varargout{3} = I;
  else
    varargout{1} = sparse(I,I,MD);      
  end  
  
  return
   
  