function varargout = assemMat_Stokes_P2P0(Mesh,EHandle,varargin)
% Assemble Galerking matrix for P2-P0 finite element discretization of Stokes problem
% \eqref{eq:VPSTaug}: piecewise quadratic continuous velocity components and piecewise
% constant pressure approximation.
%\texttt{mesh} LehrFEM mesh data structure, complete with edge information,
%Sect.~\ref{sec:meshdata} The struct MESH must at least contain the following fields:
% COORDINATES M-by-2 matrix specifying the vertices of the mesh.
% ELEMENTS N-by-3 or matrix specifying the elements of the mesh.  
% EDGES P-by-2 matrix specifying the edges of the mesh.  
% ELEMFLAG N-by-1 matrix specifying additional elementinformation. 
% \texttt{EHandle} passes function for computation of element matrix.
% See Sect.~\ref{sec:assembly} for a discussion of the generic assembly algorithm
nCoordinates = size(Mesh.Coordinates,1);
nElements = size(Mesh.Elements,1);
nEdges = size(Mesh.Edges,1);
% Preallocate memory for the efficient initialization of sparse matrix, Ex.~\ref{ex:effimpl}
I = zeros(196*nElements,1); J = zeros(196*nElements,1); A = zeros(196*nElements,1);
% Local assembly: loop over all cells of the mesh  
loc = 1:196;
for i = 1:nElements
  % Extract vertex coordinates
  vidx = Mesh.Elements(i,:);
  Vertices = Mesh.Coordinates(vidx,:);
  % Compute $14\times14$ element matrix: there are 6 local shape functions for the finite
  % element space \Blue{$\LagrFE{2}$}, and 1 (constant) local shape function for
  % \Blue{$\QFE{0}$}: \Blue{$6+6+1=13$} local shape functions for the P2-P0 scheme
  Aloc = EHandle(Vertices,Mesh.ElemFlag(i),varargin{:});
  % Add contributions to global Galerkin matrix: the numbering convention is a follows:
  % d.o.f. for \Blue{$x_1$}-components of the velocity are numbered first, then
  % \Blue{$x_2$}-components of the velocity, then the pressure d.o.f.
  eidx = [Mesh.Vert2Edge(vidx(1),vidx(2)) ...
          Mesh.Vert2Edge(vidx(2),vidx(3)) ...
          Mesh.Vert2Edge(vidx(3),vidx(1))];
  % Note: entries of an extra last row/column of the Galerkin matrix corresponding to
  % pressure d.o.f. are filled with one to enfore zero mean pressure, see
  % Ex.~\ref{ex:zeromean}
  idx = [vidx,eidx+nCoordinates,... % $\leftrightarrow$ \Blue{$v_1$}
         vidx+nCoordinates+nEdges,eidx+2*nCoordinates+nEdges, % $\leftrightarrow$ \Blue{$v_2$}
         i+2*(nEdges+nCoordinates), % $\leftrightarrow$ \Blue{$p$}
         2*(nEdges+nCoordinates)+nElements+1]; % $\leftrightarrow$ zero mean multiplier
  I(loc) = set_Rows(idx,14); J(loc) = set_Cols(idx,14); A(loc) = Aloc(:);
  loc = loc+196;
  end
  
  % Assign output arguments for creation of sparse matrix
  if(nargout > 1), varargout{1} = I; varargout{2} = J; varargout{3} = A;
  else, varargout{1} = sparse(I,J,A); end
return