function varargout = assemMat_ContrTwo2(Mesh, vHandle, varargin)
% assemMat_ContrOne Assemble topological gradient. loop over egdes
%
%   A = assemMat_ContrTwo(Mesh, vHandle, varargin) 
%   A = ASSEMMat_ContrTwo(MESH,  vHandle) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_ContrTwo(MESH, vHandle) .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemMat_ContrTwo(Mesh,vHandle);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nEdges =size(Mesh.Edges,1);
  nElements=size(Mesh.Elements,1);
  
  BndEdges=get_BdEdges(Mesh);
  
  % Preallocate memory
  
  I = zeros(2*nEdges-size(BndEdges,1),1);
  J = zeros(2*nEdges-size(BndEdges,1),1);
  A = zeros(2*nEdges-size(BndEdges,1),1);
  
  loc=1:2;
  
  %Upwinding for inner edges
  
  for i = setdiff( 1:nEdges, BndEdges)

     % Vertices
     a1=Mesh.Coordinates(Mesh.Edges(i,1),:);
     a2=Mesh.Coordinates(Mesh.Edges(i,2),:);
     
     % normal
     
     n=[a2(2)-a1(2) -a2(1)+a1(1)];
     
     % Extract nodal vectors
    
     v1 =-vHandle(a1);
     v2 =-vHandle(a2);
     
     
     % left Triangle
     
     T_l=Mesh.Edge2Elem(i,1);
     vid_l=Mesh.Elements(T_l,:);
     Vert_l=Mesh.Coordinates(vid_l,:);
     det_BK_l=abs(det([Vert_l(2,:)-Vert_l(1,:); Vert_l(3,:)-Vert_l(1,:)]));
     % sign
     p_l=-2*mod(find(vid_l==Mesh.Edges(i,2))-find(vid_l==Mesh.Edges(i,1)),3)+3;
     
     % right Triangle
     
     T_r=Mesh.Edge2Elem(i,2);
     vid_r=Mesh.Elements(T_r,:);
     Vert_r=Mesh.Coordinates(vid_r,:);
     det_BK_r=abs(det([Vert_r(2,:)-Vert_r(1,:); Vert_r(3,:)-Vert_r(1,:)]));
     % sign
     p_r=-2*mod(find(vid_r==Mesh.Edges(i,2))-find(vid_r==Mesh.Edges(i,1)),3)+3;
     % Extract global edge number
     
     % compute Matrix Entries concerning Edge i
     
     %first edge
     
     % 1.
     if ( v1*n'>=0 && v2*n'< 0 )
         B_r=p_r*abs((v1*n')^2/((v1-v2)*n'));
         B_l=p_l*abs((v2*n')^2/((v1-v2)*n'));
     end
     % 2.
     if ( v1*n'< 0 && v2*n'>= 0 )
         B_r=p_r*abs((v2*n')^2/((v1-v2)*n'));
         B_l=p_l*abs((v1*n')^2/((v1-v2)*n'));
     end
     % 3.
     if ( v1*n'>=0 && v2*n'>= 0 )
         B_r=p_r*(abs(v1*n')+abs(v2*n'));
         B_l=0;
     end
     % 4.
     if ( v1*n'< 0 && v2*n'< 0 )
         B_r=0;
         B_l=p_l*(abs(v1*n')+abs(v2*n'));
     end
     
     I(loc) = [i i];
     J(loc) = [T_l, T_r];
     A(loc) =[B_l/det_BK_l B_r/det_BK_r];
     loc = loc+2;
     
  end
  
 % Standard for boundary edges
  loc=loc(end)-1; 
 
  for i = BndEdges'
     % Vertices
     a1=Mesh.Coordinates(Mesh.Edges(i,1),:);
     a2=Mesh.Coordinates(Mesh.Edges(i,2),:);
     
     % normal
     
     n=[a2(2)-a1(2) -a2(1)+a1(1)];
     
     % Extract nodal vectors
    
     v1 =-vHandle(a1);
     v2 =-vHandle(a2);
     
     % left Triangle
     
     T_l=Mesh.Edge2Elem(i,1);
     T_r=Mesh.Edge2Elem(i,2);
     
     if (T_l==0)
         T=T_r;
         vid=Mesh.Elements(T,:);
         Vert=Mesh.Coordinates(vid,:);
         det_BK=abs(det([Vert(2,:)-Vert(1,:); Vert(3,:)-Vert(1,:)]));
         % sign
         p=-2*mod(find(vid==Mesh.Edges(i,2))-find(vid==Mesh.Edges(i,1)),3)+3;
     
         % compute Matrix Entries concerning Edge i
     
         % 1.
         if ( v1*n'>=0 && v2*n'< 0 )
          %B=p*abs((v1*n')^2/((v1-v2)*n'));
          B=p/2*(abs((v1*n')^2/((v1-v2)*n'))+abs((v2*n')^2/((v1-v2)*n')));
          %B_l=p_l*abs((v2*n')^2/((v1-v2)*n'));
         end
         % 2.
         if ( v1*n'< 0 && v2*n'>= 0 )
          %B=p*abs((v2*n')^2/((v1-v2)*n'));
          B=p/2*(abs((v1*n')^2/((v1-v2)*n'))+abs((v2*n')^2/((v1-v2)*n')));
          %B_l=p_l*abs((v1*n')^2/((v1-v2)*n'));
         end
         % 3.
         if ( v1*n'>=0 && v2*n'>= 0 )
          B=p*(abs(v1*n')+abs(v2*n'));
          %B_l=0;
         end
         % 4.
         if ( v1*n'< 0 && v2*n'< 0 )
          %B=0;
          B=p*(abs(v1*n')+abs(v2*n'));
          %B_l=p_l*(abs(v1*n')+abs(v2*n'));
         end
     else
         T=T_l;
             vid=Mesh.Elements(T,:);
         Vert=Mesh.Coordinates(vid,:);
         det_BK=abs(det([Vert(2,:)-Vert(1,:); Vert(3,:)-Vert(1,:)]));
         % sign
         p=-2*mod(find(vid==Mesh.Edges(i,2))-find(vid==Mesh.Edges(i,1)),3)+3;
     
         % compute Matrix Entries concerning Edge i
     
         % 1.
         if ( v1*n'>=0 && v2*n'< 0 )
          %B=p*abs((v1*n')^2/((v1-v2)*n'));
          B=p/2*(abs((v1*n')^2/((v1-v2)*n'))+abs((v2*n')^2/((v1-v2)*n')));
          %B=p*abs((v2*n')^2/((v1-v2)*n'));
         end
         % 2.
         if ( v1*n'< 0 && v2*n'>= 0 )
          %B=p*abs((v2*n')^2/((v1-v2)*n'));
          B=p/2*(abs((v1*n')^2/((v1-v2)*n'))+abs((v2*n')^2/((v1-v2)*n')));
          %B=p*abs((v1*n')^2/((v1-v2)*n'));
         end
         % 3.
         if ( v1*n'>=0 && v2*n'>= 0 )
          B=p*(abs(v1*n')+abs(v2*n'));
          %B=0;
         end
         % 4.
         if ( v1*n'< 0 && v2*n'< 0 )
          %B=0;
          B=p*(abs(v1*n')+abs(v2*n'));
         end     
     end
     
     I(loc) = [i];
     J(loc) = [T];
     A(loc) =[B/det_BK];
     loc = loc+1;
     
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