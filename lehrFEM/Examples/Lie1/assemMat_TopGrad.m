function varargout = assemMat_TopGrad(Mesh, varargin)
% assemMat_TopGrad Assemble topological gradient.
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
%   A = assemMat_TopGrad(Mesh,EHandle);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
 % Assign output arguments
  nE=size(Mesh.Edges,1);
  IS=zeros(nE,1);
  IZ=zeros(nE,1);
 
  for (e = 1:nE)
      IZ(e)=Mesh.Edges(e,2);
      IS(e)=Mesh.Edges(e,1);
  end
  
  A=[-ones(nE,1);ones(nE,1)];
  if(nargout > 1)
    varargout{1} = [1:nE; 1:nE]';
    varargout{2} = [IS;IZ];
    varargout{3} = A;
  else
    varargout{1} = sparse([1:nE 1:nE]',[IS;IZ],A);      
  end  
  
return  
  