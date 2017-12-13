function varargout = assemMat_Contr1f(Mesh, vHandle, varargin)
% assemMat_Contr1f assemble contraction of one-forms.
%
%   A = assemMat_Contr1f(Mesh, vHandle, varargin) 
%   A = ASSEMMat_Contr1f(MESH,  vHandle) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_Contr1f(MESH, vHandle) .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)]
%   A = assemMat_Contr1f(Mesh,V_Handle);
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
  
  I = zeros(9*nElements,1);
  J = zeros(9*nElements,1);
  A = zeros(9*nElements,1);
  
  % look for upwind nodes on each triangle
  loc=1:9;
  for i = 1:nElements
      
     B=zeros(3,3);
       
     % Vertices
     vid = Mesh.Elements(i,:);
     a1 = Mesh.Coordinates(vid(1),:);
     a2 = Mesh.Coordinates(vid(2),:);
     a3 = Mesh.Coordinates(vid(3),:);
     
     % Extract global edge numbers
    
     eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];
       
    % Determine the orientation
    
    if (Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else   p1 = -1;  end
    if (Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else   p2 = -1;  end
    if (Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else   p3 = -1;  end
     
     % Compute element mapping

     bK = a1;
     BK = [a2-bK; ...
        a3-bK];
     det_BK = abs(det(BK));
     inv_BK = inv(BK); 
  
     
     gradN1=[-1,-1]*inv_BK';
     gradN2=[1,0]*inv_BK';
     gradN3=[0,1]*inv_BK';
     
     % Extract nodal vectors
        v1=-vHandle(a1);
        v2=-vHandle(a2);
        v3=-vHandle(a3);

    % Compute decision variables Theta(T,
  
  
     % for first vertex
     y = a1 + v1;
     yhat = (y-bK)*inv_BK;
    
     if(yhat(1) >= 0 && yhat(2) >= 0)
       B(1,:)=-[0 -p2*gradN3*v1' p3*gradN2*v1'];
       if (yhat(1)==0 || yhat(2)==0)
         B(1,:)=1/2 * B(1,:);
       end   
     end
            
    % for second vertex
      
    y = a2 + v2;
    yhat = (y-bK)*inv_BK;
    
    if(yhat(2) >= 0 && sum(yhat) <= 1)
      B(2,:) = -[p1*gradN3*v2' 0 -p3*gradN1*v2'];
      if (yhat(2)==0 || sum(yhat) == 1)
        B(2,:)=1/2*B(2,:);
      end  
    end
      
    % for third index
    
    y = a3 + v3;
    yhat = (y-bK)*inv_BK;
  
    if(yhat(1) >= 0 && sum(yhat) <= 1)
      B(3,:) = -[-p1*gradN2*v3' p2*gradN1*v3' 0];
      if (yhat(1)==0 || sum(yhat) == 1)
        B(3,:)=1/2*B(3,:);
      end
    end
   
   I(loc) = set_Rows(vid,3);
   J(loc) = set_Cols(eidx,3);
   A(loc) = B(:);
   loc = loc+9;
  end
  
  
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A,nCoordinates,nEdges);      
  end
  
return