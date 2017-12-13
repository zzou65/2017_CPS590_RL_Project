function U = eval_PWDG(Mesh,omega,u,x)
%EVAL_PWDG evaluate PWDG function
%
%   UVAL = EVAL_PWDG(MESH,OMEGA,U,X) evaluates the discontinuous plane
%   wave functions given by the columns of U at the points X.
%
%   The struct MESH must contain at least the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the
%                mesh.
%    ELEMDATA     N-by-1 structure array containing at least the fields:
%       NDOFS       The number of degrees of freedom on the corresponding
%                   element.
%       DIR         P-by-2 matrix specifying the propagation directions of
%                   the plane wave basis functions.
%
%   OMEGA is the wave number of the plane wave basis functions.
%
%   U is a K-by-J matrix.  The columns of U contain vectors in the
%   space spanned by the plane wave basis functions on all the elements of
%   the mesh.
%
%   X is a I-by-2 matrix. The rows of X are coordinate pairs of the
%   evaluation points.
%
%   The function values are returned in the I-by-J matrix UVAL.  The entry
%   UVAL(i,j) contains the value of the j-th function (represented by the
%   j-th column of U) at the i-th evaluation point (ie. the i-th row of X).
%   If the i-the evaluation point lies on the boundary between multiple
%   elements, then U(i,j) contains the average of the limits of the values 
%   of the i-th function from all neighbouring elements.  If it does not
%   lie in any element, then the value is set to zero.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize
  U = zeros(size(x,1),size(u,2));
  
  % Map evaluation points to elements of the mesh
  nelem = size(Mesh.Elements,1);
  nvert = size(Mesh.Elements,2);  % Number of vertices per element
  nDofs = [Mesh.ElemData.nDofs];  % Number of degrees of freedom on each element
  nDofsSum = cumsum([0,nDofs]);   % Number of degrees of freedom on all preceading elements
  nx = size(x,1);
  
  ind = cell(nelem,1);      % list points inside each element
  elems = zeros(nx,1);      % list number of elements each point is contained in (-> averages)
  
  for j=1:nelem
    
    % Find evaluation points inside element j
    indloc = true(nx,1);
    for k=1:nvert
      dx = x - repmat(Mesh.Coordinates(Mesh.Elements(j,k),:),nx,1);
      dy = Mesh.Coordinates(Mesh.Elements(j,mod(k+1,nvert)+1),:) - Mesh.Coordinates(Mesh.Elements(j,k),:);
      edge = Mesh.Coordinates(Mesh.Elements(j,mod(k,nvert)+1),:) - Mesh.Coordinates(Mesh.Elements(j,k),:);
      normal = [-edge(2),edge(1)];
      dxn = dx(:,1)*normal(1)+dx(:,2)*normal(2);% sum(dx.*repmat(normal,nx,1),2);
      dyn = dy*normal.';
      indloc = indloc & (dxn*dyn>=0);
    end
    ind{j} = find(indloc);
    elems(ind{j}) = elems(ind{j}) + 1;
    
  end
  
  % Evaluate function
  for j=1:nelem
    if(~isempty(ind{j}))
      x0 = Mesh.Coordinates(Mesh.Elements(j,1),:);
      dir = Mesh.ElemData(j).Dir;
      loc = nDofsSum(j) + (1:nDofs(j));
      U(ind{j},:) = U(ind{j},:) + ...
        exp(i*omega*(x(ind{j},:)-x0(ones(size(ind{j})),:))*dir.')*u(loc,:);
    end
  end
  U = U./elems(:,ones(1,size(U,2)));
  
return