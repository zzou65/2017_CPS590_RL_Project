function varargout = assemMat_Contr2f(Mesh, v, varargin)
% assemMat_Contr2f assemble contraction of twoforms.
%
%   A = ASSEMMat_Contr2f(MESH,  v) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_Contr2f(MESH, v) .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)]
%   v=V_Handle(New_Mesh.Coordinates);
%   A = assemMat_Contr2f(Mesh,v);
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
  
  New_Mesh=Mesh;
  
  % Preallocate memory
  
  I = zeros(3*nElements,1);
  J = zeros(3*nElements,1);
  A = zeros(3*nElements,1);
  
  % look for upwind nodes on each triangle
  loc=1:3;
  
  rot=[0 1; -1 0];
  
  for i = 1:nElements
      
     Theta=zeros(3,2);
     Theta2=zeros(3,2);
     B=zeros(1,3);
     
     % Vertices
     vid = Mesh.Elements(i,:);
     a1 = Mesh.Coordinates(vid(1),:);
     a2 = Mesh.Coordinates(vid(2),:);
     a3 = Mesh.Coordinates(vid(3),:);
     
     % rotated egdes = normal
     
     
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
     inv_BK = inv(BK); 
  
     % Extract nodal vectors
%     
     v1 =-v(vid(1),:);
     v2 =-v(vid(2),:);
     v3 =-v(vid(3),:);
%      v1 =-v(a1);
%      v2 =-v(a2);
% %      v3 =-v(a3);
%  
%      % endpoints of nodal vectors
     y1 = a1 + v1;
     y2 = a2 + v2;
     y3 = a3 + v3;
%      
%      % and their counterparts on the reference triangle
%      
%      yhat1 = (y1-bK)*inv_BK;
%      yhat2 = (y2-bK)*inv_BK;
%      yhat3 = (y3-bK)*inv_BK;
%       
%       yhat1=yhat1/norm(yhat1);
%       yhat2=yhat2/norm(yhat2);
%       yhat3=yhat3/norm(yhat3);
% %      % Compute decision variables Theta
%      
%      % for first edge
%      if ((yhat2-[1 0])*[1;1]<-TOL*sqrt(2))
%          Theta(1,1)=1;
%      end
%      if ((yhat3-[0 1])*[1;1]<-TOL*sqrt(2))
%          Theta(1,2)=1;
%      end
%      
%      % for second edge
%      if ((yhat3-[0 1])*[-1;0]<-TOL)
%          Theta(2,1)=1;
%      end
%      if ((yhat1)*[-1;0]<-TOL)
%          Theta(2,2)=1;
%      end
%     
%      % for third edge
%      if ((yhat1)*[0;-1]<-TOL)
%          Theta(3,1)=1;
%      end
%      if ((yhat2-[1 0])*[0;-1]<-TOL)
%          Theta(3,2)=1;
%      end

     % Compute decision variables Theta
     
     % for first edge
     if ((v2)*n1'>TOL*norm(n1)*norm(v2))
         Theta(1,1)=1;
     end
     if ((v3)*n1'>TOL*norm(n1)*norm(v3))
         Theta(1,2)=1;
     end
     
     % for second edge
     if ((v3)*n2'>TOL*norm(v3)*norm(n2))
         Theta(2,1)=1;
     end
     if ((v1)*n2'>TOL*norm(v1)*norm(n2))
         Theta(2,2)=1;
     end
    
     % for third edge
     if ((v1)*n3'>TOL*norm(v1)*norm(n3))
         Theta(3,1)=1;
     end
     if ((v2)*n3'>TOL*norm(v2)*norm(n3))
         Theta(3,2)=1;
     end
     
     % compute Matrix Entries concerning Element i
     
     %first edge
     
     % 1.
     if ( Theta(1,1)==1 && Theta(1,2)==0 )
         %B(1,1)=p1*abs((v2*n1')^2/((v2-v3)*n1'));
         B(1,1)=p1*abs(v2*n1');
     end
     % 2.
     if (Theta(1,1)==0 && Theta(1,2)==1)
         %B(1,1)=p1*abs((v3*n1')^2/((v2-v3)*n1'));
         B(1,1)=p1*abs(v3*n1')
     end
     % 3.
     if (Theta(1,1)==1 && Theta(1,2)==1)
         B(1,1)=p1*(abs(v2*n1')+abs(v3*n1'));
     end
     % 4.
     % if (Theta(1,1)=0 && Theta(1,2)=0)
     % B(1,1)=
     % end
     
     
     % second edge
     
     % 1.
     if (Theta(2,1)==1 && Theta(2,2)==0)
         B(1,2)=p2*abs((v3*n2')^2/((v1-v3)*n2'));
         %B(1,2)=p2*abs(v3*n2');
     end
     % 2.
     if (Theta(2,1)==0 && Theta(2,2)==1)
         B(1,2)=p2*abs((v1*n2')^2/((v1-v3)*n2'));
         %B(1,2)=p2*abs(v1*n2')^2;
     end
     % 3.
     if (Theta(2,1)==1 && Theta(2,2)==1)
         B(1,2)=p2*(abs(v1*n2')+abs(v3*n2'));
     end
     % 4.
     % if (Theta(2,1)=0 && Theta(2,2)=0)
     % B(1,2)=
     % end
     
     % Third edge
     
     % 1.
     if (Theta(3,1)==1 && Theta(3,2)==0)
         B(1,3)=p3*abs((v1*n3')^2/((v1-v2)*n3'));
         %B(1,3)=p3*abs(v1*n3');
     end
     % 2.
     if (Theta(3,1)==0 && Theta(3,2)==1)
         B(1,3)=p3*abs((v2*n3')^2/((v1-v2)*n3'));
         %B(1,3)=p3*abs(v1*n3');
     end
     % 3.
     if (Theta(3,1)==1 && Theta(3,2)==1)
         B(1,3)=p3*(abs(v2*n3')+abs(v1*n3'));
     end
     % 4.
     % if (Theta(3,1)=0 && Theta(3,2)=0)
     % B(1,3)=
     % end

     I(loc) = eidx;
     J(loc) = [i,i,i];
     A(loc) =B(:)./(det_BK);
     loc = loc+3;
     
  end
  
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A,nEdges,nElements);      
  end
  
return