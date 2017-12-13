function varargout = assemMat_Contr1f_QFE(Mesh, vHandle, varargin)
% assemMat_Contr1f_QFE assemble contraction of.
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
  
  I = zeros(18*nElements,1);
  J = zeros(18*nElements,1);
  A = zeros(18*nElements,1);
  
  % look for upwind nodes on each triangle
  loc=1:18;
  for i = 1:nElements
      
     B=zeros(6,3);
       
     % Vertices
     vid = Mesh.Elements(i,:);
     a1 = Mesh.Coordinates(vid(1),:);
     a2 = Mesh.Coordinates(vid(2),:);
     a3 = Mesh.Coordinates(vid(3),:);
     
     a4=(a2+a3)/2;
     a5=(a3+a1)/2;
     a6=(a1+a2)/2;
     
     % Extract global edge numbers
    eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];
     idx=eidx+nCoordinates;   
       
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
        v4=-vHandle(a4);
        v5=-vHandle(a5);
        v6=-vHandle(a6);

    % Compute element entries according to upwind definition
  
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
      B(3,:) =-[-p1*gradN2*v3' p2*gradN1*v3' 0];
      if (yhat(1)==0 || sum(yhat) == 1)
        B(3,:)=1/2*B(3,:);
      end
    end
    
  % first midpoint     
  yhat = (a4+v4-bK)*inv_BK;
  if(sum(yhat) <= 1)
    B(4,:) = -[p1/2*(gradN3-gradN2)*v4' p2/2*gradN1*v4' -p3/2*gradN1*v4'];
    end
  if (sum(yhat)==1)
     B(4,:)=1/2* B(4,:);    %two elements contribute to matrix 
  end
  
  % second midpoint    
  yhat = (a5+v5-bK)*inv_BK;
  if(yhat(1) >= 0)
    B(5,:) =-[-p1/2*gradN2*v5' p2/2*(gradN1-gradN3)*v5' p3/2*gradN2*v5']; 
  end
  if (yhat(1)==0)
    B(5,:)=1/2* B(5,:);    %two elements contribute to matrix 
  end

  % third midpoint 
  yhat = (a6+v6-bK)*inv_BK;
  if(yhat(2) >= 0)
    B(6,:) =-[p1/2*gradN3*v6'  -p2/2*gradN3*v6' p3/2*(gradN2-gradN1)*v6'];
  end
  if (yhat(2)==0)
     B(6,:)=1/2*B(6,:);    %two elements contribute to matrix 
  end
   I(loc) = set_Rows([vid idx],3);
   J(loc) = set_Cols(eidx ,6);
   A(loc) = B(:);
   loc = loc+18;
  end
  
  
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A,nCoordinates+nEdges,nEdges);      
  end
  
return
