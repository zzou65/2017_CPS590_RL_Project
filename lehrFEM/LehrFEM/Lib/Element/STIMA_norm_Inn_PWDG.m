function Aloc = STIMA_norm_Inn_PWDG(Edge,Normal,EdgeData,LData,RData,varargin)
%STIMA_NORM_INN_PWDG Inner edge contributions to DG norm for plane waves
%   
%   ALOC = STIMA_NORM_INN_PWDG(EDGE,NORMAL,EDGEDATA,LDATA,RDATA)
%   computes the contributions to the DG-norm inner product from inner
%   edges for plane waves.
%  
%   EDGE is 2-by-2 matrix whose rows contain the start and end nodes of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the unit normal with respect to
%   the current edge EDGE.  The orientation of the normal vector should be
%   from the right element to the left element.
%
%   EDGEDATA is a structure that contains at least the fields:
%    NDOFS    Total number of degrees of freedom on elements adjacent to
%             current edge.
%    L2       The L2 inner product matrix on the current edge.
%    LENGTH   Length of edge.
%
%   The structs LDATA and RDATA conatin the left and right hand side
%   element data:
%    ELEMENT  Integer specifying the neighbouring element.
%    ELEMDATA Structure containing at least the fields:
%       NDOFS   The number of degrees of freedom on the current element.
%       DIR     A P-by-2 matrix containing the propagation directions of
%               the plane wave basis functions in its rows.
%    VERTICES 3-by-2 or 4-by-2 matrix specifying the vertices of the
%             neighbouring element.
%    EDGELOC  Integer specifying the local edge number on the neighbouring
%             element.
%    MATCH    Integer specifying the relative orientation of the edge with
%             respect to the orientation of the neighbouring element.
%
%   See also STIMA_norm_Bnd_PWDG, STIMA_Lapl_Vol_PWDG, MASS_Vol_PWDG,
%   assemMat_Inn_PDG

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize
  nDofsL = LData.ElemData.nDofs;    % Number of basis functions on left element
  nDofsR = RData.ElemData.nDofs;    % Number of basis functions on right element
  nDofs = EdgeData.nDofs;           % Total number of basis functions on left and right elements
  h = EdgeData.Length;              % Length of edge
  
  % Values will be represented as linear combinations of plane waves
  % with direction vectors Dir and value one at Edge(1,:)
  
  % Calculate traces from the left
  Lv = ones(nDofsL,1);
  
  % Calculate traces from the right
  Rv = ones(nDofsR,1);
  
  % Compute jump
  vN = [-Lv;Rv];
  
  % The sought integral is the pointwise product of a factor and the L2
  % inner product matrix on the current edge.
  
  % Calculate integrand at first endpoint of edge
  Aloc = h*transpose(vN(:,ones(nDofs,1))).*conj(vN(:,ones(nDofs,1)));
  
  % Multiply by L2 inner product matrix (mass matrix)
  Aloc = Aloc.*EdgeData.L2;
  
return