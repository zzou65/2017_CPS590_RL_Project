function varargout = set_Data_PWDG(Mesh,varargin)
%SET_DATA_PWDG Set parameters for plane wave discontinuous Galerkin method
%
%   MESH = SET_DATA_PWDG(MESH,OPTIONS) creates or updates information on
%   the discontinuous plane wave discretization in the mesh data structure
%   MESH.  The options OPTIONS may be given either as a data structure or
%   as a list of field names and values, see below for details.  In the
%   latter case, only the fields whose values differ from the default or
%   previous values need to be specified.
%
%   [MESH,DATA] = SET_DATA_PWDG(MESH,OPTIONS) additionally returns a struct
%   DATA containg information on the discretization.
%
%   OPTIONS = SET_DATA_PWDG with no arguments returns the default options.
%
%   The option parameters are:
%
%   Elements: list of elements to update, values:
%     [] [default]: update all elements
%     (vector of length n): only the elements listed in the array are
%         updated.  Other information is required to have the same array
%         form where indicated.
%
%   Dir: normalized propagation direction vectors for plane wave basis
%       functions, values:
%     (p-by-2 matrix): each row represents one propagation direction; all
%         are used on every element.
%     (n-by-2 matrix): each row represents one propagation direction.  The
%         i-th row is used on element Elements(i).  This form is required
%         if Elements is nonempty.
%
%   nDofs: number of new degrees of freedom per element, values:
%     (scalar): same value is used on each element.
%     (vector): if Elements is empty, the i-th term of nDofs is used on the
%         i-th element listed in the mesh data structure, otherwise it
%         corresponds to the i-th element in the Elements argument.
%
%   RemDir: degrees of freedom that are to be removed, values:
%     (vector of length n): if Elements is empty, the n degrees of freedom
%         listed are removed from every element.  Otherwise, the n-th term
%         is removed on the n-th element in Elements.
%
%   DofData: arbitrary data corresponding to degrees of freedom, values:
%     (row vector): this vector is stored for all new degrees of freedom.
%     (m-by-k matrix): the i-th row is used for the i-th new degree of
%         freedom, counted locally on each element if Elements is empty and
%         globally of it isn't.
%
%   MakeDir: method to use for generating new propagation directions,
%       values:
%     'uniform' [default]: new propagation directions are distributed
%         equidistantly on the unit circle, with [1 0] as the first vector.
%     'orient': new propagation directions are distributed equidistantly on
%         the unit circle, with the first vector determined by the angle
%         Phi0.
%     'rand0': new propagation directions are distributed equidistantly
%         with a random first vector.
%     'random': new propagation directions are chosen randomly.
%
%   Phi0: angle between first propagation direction and [1 0], values:
%     (scalar): this value is used on all elements.
%     (vector): the i-th value is used on the i-th element of the mesh if
%         Elements is empty and on the i-th element listed in Elements if
%         the latter is not empty.
%
%   NewDir: add or replace propagation directions, values:
%     'add' [default]: append new propagation directions to existing
%         values, if values are already stored in the mesh data structure.
%     'replace': overwrite old values.  The option Elements must be set to
%         [] in this case, ie. degrees of freedom are replaced on all
%         elements.
%
%   OrderDirs: order direction vectors, values:
%     true [default]: order propagation directions counterclockwise along 
%         the unit circle, starting in the x-direction.
%     false: store propagation directions in arbitrary order.
%
%   Omega: wave number of Helmholtz equation, values:
%     (arbitrary): passed as the first argument to all flux parameters
%         given as function handles.
%     Note that a scalar value for Omega must be given in order for the L2
%     inner product matrices on edges to be precomputed.  These matrices
%     are used in various assembly routines.
%
%   FluxParams: list of names of flux parameters, values:
%     (cell array of strings) [default: {'a','b','c','d'}].
%     Note that, even if not all flux parameters are changed, all of them
%     must be contained in this list.  Parameters that are removed from
%     this option are deleted.
%
%   a: scalar coefficient for the (jump-jump) penalty term in the
%       numerical flux for the function value, values:
%     (scalar) [default: 0.5]: use value on all edges.
%     (function handle): use the value returned by the function; it should
%         take the arguments Omega, h, p, where Omega is the option
%         discussed above, h is the edge length and p is the maximal number
%         of plane waves on a neighbouring element.
%
%   b: scalar coefficient for a term containing jumps of the normal 
%       derivative in the numerical flux for the gradient, values:
%     (scalar) [default: 0.5]: use value on all edges.
%     (function handle): use the value returned by the function; it should
%         take the arguments Omega, h, p, where Omega is the option
%         discussed above, h is the edge length and p is the maximal number
%         of plane waves on a neighbouring element.
%
%   c: vector flux coefficient, values:
%     (1-by-2 vector) [default: [0 0]]: use value on all edges.
%     (function handle): use the value returned by the function; it should
%         take the arguments Omega, h, p, where Omega is the option
%         discussed above, h is the edge length and p is the maximal number
%         of plane waves on a neighbouring element.
%
%   d: scalar coefficient for flux on impedance boundary, should be between
%       0 and 1, values:
%     (scalar) [default: 0.5]: use value on all edges.
%     (function handle): use the value returned by the function; it should
%         take the arguments Omega, h, p, where Omega is the option
%         discussed above, h is the edge length and p is the maximal number
%         of plane waves on a neighbouring element.
%
%   Post: function handles for postcomputations on mesh, values:
%     (function handle): function is applied to mesh after all the data is
%         added to the mesh data structure.  The function should take as
%         argument the mesh data structure and return a (modified) mesh
%         data structure.
%     (cell array): each element of the cell array should be a function
%         handle; they are applied to the mesh after all the data is added
%         to the mesh data structure in the order in which they are listed.
%     Note that errors can occure if certain parts of the mesh data
%     structure are modified, such as the vertex coordinates or precomputed
%     inner products.  Flux parameters may be changed without limitations.
%
%   The input argument MESH must contain at least the following fields:
%
%   Coordinates: M-by-2 matrix specifying the vertices of the mesh.
%
%   Elements: N-by-3 or N-by-4 matrix specifying the elements of the 
%       mesh.
%
%   Edges: L-by-2 matrix specifying all edges of the mesh.
%
%   Edge2Elem: L-by-2 matrix connecting edges to elements. The first column
%       specifies the left hand side element and the second column
%       specifies the right hand side element.
%
%
%   Three struct arrays are added to the mesh data structure MESH,
%   containing information relevant to the full mesh, the elements of the
%   mesh and the edges of the mesh.
%
%   MESH.Data stores parameters relevant to the whole mesh.  These are:
%
%   Options: a structure containg all of the fields listed as options above
%       with values corresponding to the mesh data structure MESH.  This
%       structure, with modified values, may be used as an argument in this
%       function.  Note that the fields Elements, RemDir, Dir and nDofs are
%       reset to their default value [].
%
%   MESH.ElemData(i) stores data specific to the i-th element of the mesh.
%   Its fields are:
%
%   Area: the area of the element.
%
%   Diameter: the diameter of the element.
%
%   nDofs: number of basis functions on the current element.
%
%   Dir: nDofs-by-2 matrix containg the normalized propagation directions 
%       of the plane wave basis functions in its rows.
%
%   Phi: nDofs-by-1 vector containg the angles (between 0 and 2*PI) of the
%       plane wave propagation directions.
%
%   DofData: nDofs-by-k matrix containing arbitrary data corresponding to
%       the individual degrees of freedom.
%
%   Ind: 1-by-nDofs vector containing the global indices of the degrees of
%       freedom on the current element
%
%   MESH.EdgeData(j) stores parameters for the j-th edge of the mesh.  It
%   has the following fields:
%
%   Length: length of the edge.
%
%   nDofs: total number of degrees of freedom on all neighbouring elements.
%
%   sameDofs: logical specifying whether or not the degrees of freedom on
%       the neighbouring elements are identical.
%
%   (param) where param is an element of the field OPTIONS.FluxParams (see
%       above): numerical value of the flux parameter param on the current
%       edge.
%
%   L2: the L2 inner product matrix on the current edge.  This is used in
%       various assembly routines.
%
%
%   The optional second output argument DATA is a structure containing the
%   following fields:
%
%   Options: a structure containg all of the fields listed as options above
%       with values corresponding to the mesh data structure MESH.  This
%       structure, with modified values, may be used as an argument in this
%       function.  Note that the fields Elements, RemDir, Dir and nDofs are
%       reset to their default value [].
%
%   I0, I1: indices for basis transformation.  If U0 is a vector in the
%       original basis, before applying this function to the mesh, and U1
%       represents the same vector in the new basis, up to the basis
%       functions that were removed, then the indices of nonzero entries of
%       U1 are I1 and U1(I1) = U0(I0).
%
%
%   Examples:
%
%   1) Create uniform basis functions
%
%     a. To initialize or add five equidistantly distributed basis
%     functions to each element,
%
%       Mesh = set_Data_PWDG(Mesh,'nDofs',5,'Omega',pi);
%
%     b. To set the basis functions on all elements to plane waves with
%     propagation directions in the four coordinate directions and south
%     west,
%
%       d = [-1 -1]/sqrt(2);
%       Dir = [1 0; 0 1; -1 0; 0 -1; d];
%       Mesh = set_Data_PWDG(Mesh,'Dir',Dir,'NewDir','replace','Omega',pi);
%
%   2) Selectively modify basis functions
%
%     a. To add the plane wave with propagation direction Dir(i,:) to the
%     element Elements(i),
%
%       Mesh = set_Data_PWDG(Mesh,'Elements',Elements,'Dir',Dir);
%
%     b. To add the plane wave with propagation direction DirNew(i,:) to
%     the element ElemNew(i) and remove the DirRem(j)-th basis function
%     from the element ElemRem(j),
%
%       Mesh = set_Data_PWDG(Mesh,'Elements',ElemRem,'RemDir',DirRem);
%       Mesh = set_Data_PWDG(Mesh,'Elements',DirNew,'Dir',DirNew);
%
%     c. To add basis functions as in 2)a. and transform a vector U to the
%     new basis,
%
%       [Mesh,Data] = set_Data_PWDG(Mesh,'Elements',Elements,'Dir',Dir);
%       U_old = U;
%       U = zeros(size(Mesh.Elements,1),1);
%       U(Data.I1) = U_old(Data.I0);
%
%   3) Set flux parameters
%     
%     a. To use constant values for the flux parameters,
%
%       Mesh = set_Data_PWDG(Mesh,'a',1,'b',0.2,'c',[1 0],'d',0.75);
%
%     b. To use mesh-dependant values for a and b, and default values for c
%     and d,
%
%       Mesh = set_Data_PWDG(Mesh,...
%         'a',@(omega,h,varargin) 2/(omega*h),...
%         'b',@(omega,h,varargin) 0.1*omega*h,...
%         'Omega',12*pi);
%
%   4) Alternative syntax
%
%     a. Create default option structure and modify desired values
%
%       Options = set_Data_PWDG();
%       Options.nDofs = 7;
%       Mesh = set_Data_PWDG(Mesh,Options);
%
%     b. Use option structure from previous application of function and
%     modify desired values
%
%       Options = Mesh.Data.Options;
%       Options.nDofs = 1;
%       Options.MakeDir = 'random';
%       Mesh = set_Data_PWDG(Mesh,Options);
%

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%% Options

  % Return default options if no arguments are given
  if(nargin==0)

    % Discretization arguments
    Options.Elements = []; % or vector of relevant elements
    Options.Dir = []; % p-by-2 matrix or n-by-2 matrix with n = size(Elements,1)
    Options.nDofs = []; % scalar or vector of length = number of elements
    Options.RemDir = []; % vector of length n, corresponding to Elements
    Options.DofData = []; % Data corresponding to new degrees of freedom
    Options.MakeDir = 'uniform'; % or 'rand0' or 'random' or 'orient'
    Options.Phi0 = 0; % or some number
    Options.NewDir = 'add'; % or 'replace'
    Options.OrderDirs = true; % order direction vectors along S1
    
    Options.Omega = []; % or scalar
    
    % Flux parameters
    Options.FluxParams = {'a','b','c','d'}; % list of flux parameters
    Options.a = 0.5; % or a(omega,h,p,...)
    Options.b = 0.5; % or b(omega,h,p,...)
    Options.c = [0 0]; % or c(omega,h,p,...)
    Options.d = 0.5; % or d(omega,h,p,...)
    
    % Other parameters
    Options.Post = {}; % function handles for postcomputations on mesh
    
    % Assign options as output argument
    varargout{1} = Options;
    
    return
  end
  
  % Read options
  if(length(varargin)==1)
    Options = varargin{1};
  else
    if(isfield(Mesh,'Data') && isfield(Mesh.Data,'Options'))
      Options = Mesh.Data.Options;
    else
      Options = set_Data_PWDG;
    end
    for i=1:floor(0.5*length(varargin))
      Options.(varargin{2*i-1}) = varargin{2*i};
    end
  end


%% Initialize
  
  nElements = size(Mesh.Elements,1);% number of elements in mesh
  nEdges = size(Mesh.Edges,1);      % number of edges in mesh
  nVert = size(Mesh.Elements,2);    % number of vertices / edges per element
  
  % Initialize data structures
  if(~isfield(Mesh,'EdgeData') || ~isequal(size(Mesh.EdgeData,1),nEdges))
    Mesh.EdgeData = repmat(struct(),nEdges,1);
  end
  if(~isfield(Mesh,'ElemData') || ~isequal(size(Mesh.ElemData,1),nElements))
    Mesh.ElemData = repmat(struct(),nElements,1);
  end
  if(~isfield(Mesh,'Data'))
    Mesh.Data = struct();
  end
  
  % Initialize second output argument and temporary data structure
  Data = struct;
  LocData = repmat(struct(),nElements,1);

  % Initialize basis transformation:
  %   If u0 is vector on original element and u1 its representation in
  %   the modified basis, then u0(I0) = u1(I1), nDofs0 is the number of
  %   original basis fn, nDofs1 is the size of the new basis and nInd is
  %   the size of the intersection of the to basis sets.
  if(isfield(Mesh.ElemData,'nDofs'))
    for i=1:nElements
      LocData(i).I0 = 1:Mesh.ElemData(i).nDofs;
      LocData(i).I1 = 1:Mesh.ElemData(i).nDofs;
      LocData(i).nDofs0 = Mesh.ElemData(i).nDofs;
      LocData(i).nDofs1 = LocData(i).nDofs0;
      LocData(i).nInd = LocData(i).nDofs0;
    end
  else
    [LocData.I0] = deal(zeros(1,0));
    [LocData.I1] = deal(zeros(1,0));
    [LocData.nDofs0] = deal(0);
    [LocData.nDofs1] = deal(0);
    [LocData.nInd] = deal(0);
  end
  [LocData.NewInd] = deal(zeros(1,0));      % indices of new basis functions
  [LocData.nNewDofs] = deal(0);
  [LocData.updateL2] = deal(false);
  
  % Compute geometrical data for mesh
  if(~isfield(Mesh.EdgeData,'Length'))
    for i=1:nEdges
      % Compute edge length
      Edge = Mesh.Coordinates(Mesh.Edges(i,:),:);
      Mesh.EdgeData(i).Length = norm(Edge(2,:)-Edge(1,:));
    end
  end
  
  if(~isfield(Mesh.ElemData,'Area'))
    for i=1:nElements
      
      % Get coordinates of vertices
      Vertices = Mesh.Coordinates(Mesh.Elements(i,:),:);
      x = Vertices(:,1);
      y = Vertices(:,2);
      
      % Calculate area of element (=diagonal entries) using formula due to
      % Meister & Gauss (cf. for example Wikipedia entry on polygons).
      dy = y(mod(1:nVert,nVert)+1) - y(mod(-1:nVert-2,nVert)+1);
      Mesh.ElemData(i).Area = 0.5*sum(x.*dy);
      
      % Compute diameter of element
      upper = triu(true(nVert),1);
      [X1,X2] = meshgrid(x,x);
      X = X1(upper) - X2(upper);
      [Y1,Y2] = meshgrid(y,y);
      Y = Y1(upper) - Y2(upper);
      D = X.^2 + Y.^2;
      Mesh.ElemData(i).Diameter = sqrt(max(D));      
      
    end
  end
  
%% Handle elements
  % Update degrees of freedom
  
  % Initialize
  
  % Find out what to do
  doRemove = ~isempty(Options.RemDir);
  doCreateDir = isempty(Options.Dir) && ~isempty(Options.nDofs);
  doAllElem = isempty(Options.Elements);
  doNew = doCreateDir || ~isempty(Options.Dir);
  
  % Find indices of edges that need to be modified
  if(doAllElem)
    ElementSet = 1:nElements;
  else % indices of equal entries in Element can be constructed from K and numElem.
    [ElementSort,K] = sort(Options.Elements);
    [ElementSet,I] = unique(ElementSort);
    numElem = diff([0;I(:)]);
  end
  
  % Make nDofs into a vector, to simplify code
  if(isscalar(Options.nDofs))
    nDofs = Options.nDofs(ones(nElements,1));
  else
    nDofs = Options.nDofs;
  end
  
  % Change size of DofData to simplify code
  if(doNew && isempty(Options.DofData))     % no data is given
    if(doAllElem)
      Options.DofData = zeros(1,0);
    elseif(doCreateDir)
      Options.DofData = zeros(sum(nDofs),0);
    else
      Options.DofData = zeros(size(Options.Dir,1),0);
    end
  end
  if(doNew && size(Options.DofData,1)==1)   % a single set of data is given for all new dofs
    if(doCreateDir && doAllElem)          
      Options.DofData = Options.DofData(ones(Options.nDofs,1),:);
    elseif(doCreateDir && ~doAllElem)
      Options.DofData = Options.DofData(ones(sum(nDofs),1),:);
    else
      Options.DofData = Options.DofData(ones(size(Options.Dir,1),1),:);
    end
  end
  
  % Reset data if necessary
  doClear = ~isfield(Mesh.ElemData,'Dir') ...
      || (isequal(Options.NewDir,'replace') && doNew);
  if(doClear)
    [Mesh.ElemData.Dir] = deal(zeros(0,2));
    [Mesh.ElemData.Phi] = deal(zeros(0,1));
    [Mesh.ElemData.DofData] = deal(zeros(0,size(Options.DofData,2)));
    [Mesh.ElemData.nDofs] = deal(0);
    [LocData.nDofs1] = deal(0);
  end
  
  % Check if omega has changed or if it is not given at all
  doNewOmega = doClear;
  if(isfield(Mesh.Data,'Options'))
    doNewOmega = doNewOmega || ~isequal(Options.Omega,Mesh.Data.Options.Omega);
  end
  skipOmega = isempty(Options.Omega);
  if(doNewOmega || skipOmega)
    [Mesh.EdgeData.L2] = deal([]);
    
    % Update basis transformation and other local data
    [LocData.I0] = deal(zeros(1,0));
    [LocData.I1] = deal(zeros(1,0));
    [LocData.nInd] = deal(0);
    [LocData.updateL2] = deal(true);
  end
  
  
  % Define routine for creating new direction vectors
  if(doCreateDir)
    switch Options.MakeDir
      case 'uniform'
        makePhi = @(p,varargin) 2*pi*(0:1/p:1-1/p)';
      case 'orient'
        if(isa(Options.Phi0,'numeric'))
          makePhi = @(p,varargin) Options.Phi0(min(varargin{1},length(Options.Phi0))) ...
            + 2*pi*(0:1/p:1-1/p)';
        elseif(isa(Options.Phi0,'function_handle'))
          makePhi = @(p,varargin) Options.Phi0(p,varargin{:}) ...
            + 2*pi*(0:1/p:1-1/p)';
        end
      case 'rand0'
        makePhi = @(p,varargin) 2*pi*(rand(1)+(0:1/p:1-1/p)');
      case 'random'
        makePhi = @(p,varargin) 2*pi*rand(p,1);
    end
    getCoordinates = @(phi) [cos(phi),sin(phi)];
    makeDir = @(p,varargin) getCoordinates(makePhi(p,varargin{:}));
  end
  
  % Create direction vectors
  j0 = 0;
  for i0=1:length(ElementSet)
    
    i = ElementSet(i0);
    
    % Extract current direction vectors
    Dir = Mesh.ElemData(i).Dir;
    Phi = Mesh.ElemData(i).Phi;
    DofData = Mesh.ElemData(i).DofData;
    
    % Find indices of current element in element list
    if(~doAllElem)
      j1 = j0 + numElem(i0);
      k = K(j0+1:j1);
      j0 = j1;
    end
    
    % Remove degrees of freedom on current element
    if(doRemove)
      if(doAllElem) % same thing on each element
        keep = setdiff(1:Mesh.ElemData(i).nDofs,Options.RemDir);
      else % different dofs on each element
        keep = setdiff(1:Mesh.ElemData(i).nDofs,Options.RemDir(k));
      end
      Dir = Dir(keep,:);
      Phi = Phi(keep);
      DofData = DofData(keep,:);
      
      % Update basis transformation
      LocData(i).I0 = LocData(i).I0(keep);
      LocData(i).nDofs1 = length(keep);
      LocData(i).I1 = 1:LocData(i).nDofs1;
      LocData(i).nInd = LocData(i).nDofs1;
    end
    
    % Add degrees of freedom
    if(doNew)
      % Construct new degrees of freedom on current element
      if(doCreateDir && doAllElem) % create direction vectors
        NewDir = makeDir(nDofs(i),i);
        NewDofData = Options.DofData;
      elseif(doCreateDir) % create direction vectors
        NewDir = makeDir(nDofs(k),k);
        NewDofData = Options.DofData(k,:);
      elseif(doAllElem) % use same direction vectors for all elements
        NewDir = Options.Dir;
        NewDofData = Options.DofData;
      else % extract relevant direction vectors from list
        NewDir = Options.Dir(k,:);
        NewDofData = Options.DofData(k,:);
      end
      
      % Append direction vectors
      Dir = [Dir;NewDir];   % M-Lint warning is irrelevant
      
      % Compute angles for new direction vectors
      NewPhi = getPhi(NewDir);
      Phi = [Phi;NewPhi];   % M-Lint warning is irrelevant
      
      % Append data for new degrees of freedom
      DofData = [DofData;NewDofData]; % M-Lint warning is irrelevant
      
      % Update basis transformation
      LocData(i).nDofs1 = size(Dir,1);
      LocData(i).nNewDofs = size(NewDir,1);
      LocData(i).NewInd = LocData(i).nDofs1-LocData(i).nNewDofs+1:LocData(i).nDofs1;
    end
    
    % Order direction vectors
    if(Options.OrderDirs)
      [Phi,ind] = sort(Phi);
      Dir = Dir(ind,:);
      DofData = DofData(ind,:);
      
      % Update basis transformation
      [nums,ind1] = sort(ind);
      LocData(i).I1 = ind1(LocData(i).I1);
      LocData(i).NewInd = ind1(LocData(i).NewInd);
    end
    
    % Update data structure
    Mesh.ElemData(i).nDofs = size(Dir,1);
    Mesh.ElemData(i).Dir = Dir;
    Mesh.ElemData(i).Phi = Phi;
    Mesh.ElemData(i).DofData = DofData;
    
  end
  
  % Constuct global indices for local basis functions
  if(doNew || doRemove)
    nDofSum = cumsum([0,Mesh.ElemData.nDofs]);
    for i=1:nElements
      Mesh.ElemData(i).Ind = (nDofSum(i)+1):nDofSum(i+1);
    end
    
    % Set flags to update L2 inner product
    [LocData(ElementSet).updateL2] = deal(true);
  end
  
  % If omega was changed, mark all indices as new
  if(doNewOmega)
    for i=1:nElements
      LocData(i).NewInd = 1:Mesh.ElemData(i).nDofs;
      LocData(i).nNewDofs = Mesh.ElemData(i).nDofs;
    end
  end
  

%% Handle edges
  % Compute flux parameters and edge data
  
  FluxParams = Options.FluxParams;
  
  % Distinguish between flux parameters given as function handles and those
  % given as constants
  fhandle_flag = cellfun(@(nm)isa(Options.(nm),'function_handle'),FluxParams);
  const_flag = ~fhandle_flag;
  fhandles = find(fhandle_flag);
  consts = find(const_flag);
  
  % Check if parameters have been changed
  doFluxes = ~isfield(Mesh.Data,'Options');
  if(~doFluxes)
    if(~isequal(FluxParams,Mesh.Data.Options.FluxParams))
      remOpts = setdiff(Mesh.Data.Options.FluxParams,FluxParams);
      Options = rmfield(Options,remOpts);
      for j=1:length(remOpts)
        if(isfield(Mesh.EdgeData,remOpts{j}))
          Mesh.EdgeData = rmfield(Mesh.EdgeData,remOpts{j});
        end
      end
      doFluxes = true;
    else
      for j=1:length(FluxParams)
        if(~isequal(Options.(FluxParams{j}),Mesh.Data.Options.(FluxParams{j})))
          doFluxes = true;
          break;
        end
      end
    end
  end
  
  % Set constant flux values
  if(doFluxes && any(const_flag))
    for j=1:length(consts)
      param = FluxParams{consts(j)};
      [Mesh.EdgeData.(param)] = deal(Options.(param));
    end
  end
  
  % Find edges that need updating
  if(doAllElem || doFluxes || doNewOmega)
    EdgeSet = 1:nEdges;
  else
    ElementSetVert = Mesh.Elements(ElementSet,:);
    ind = [nVert,1:nVert-1];
    Edges = Mesh.Vert2Edge(ElementSetVert,ElementSetVert(:,ind));
    Edges = unique(Edges)';
    EdgeSet = Edges(Edges~=0);
  end
  
  % Update data on edges
  doFunHandles = any(fhandle_flag);
  for i=EdgeSet
    
    % Update edge information on degrees of freedom of neighbouring
    % elements
    ind = Mesh.Edge2Elem(i,:);
    ind = ind(ind~=0);              % indices of elements adjacent to edge
    p0 = [Mesh.ElemData(ind).nDofs];
    nDofs = sum(p0);
    Mesh.EdgeData(i).nDofs = nDofs;
    if(length(ind)==2)
      Mesh.EdgeData(i).sameDofs = ...
        isequal(Mesh.ElemData(ind(1)).Dir,Mesh.ElemData(ind(2)).Dir);
    end
    
    % Evaluate function handles for fluxes
    if(doFunHandles)
      omega = Options.Omega;
%       h = Mesh.EdgeData(i).Length;
      h = min([Mesh.ElemData(ind).Diameter]);
      p = max(p0);
      for j=1:length(fhandles)
        param = FluxParams{fhandles(j)};
        Mesh.EdgeData(i).(param) = Options.(param)(omega,h,p);
      end
    end
    
%% Precompute L2 inner product on edge
    
    % Check if anything changed
    if(~any([LocData(ind).updateL2]) || skipOmega)
      continue
    end
    
    % The product of two plane waves is again a plane wave, with the
    % difference (because of complex conjugation) of the propagation
    % directions as propagation direction.  Its integral is composed of the
    % following terms:
    %   expX0: the value at the first endpoint of the edge.
    %   Length: the length of the edge.
    %   expBmA: (exp(z)-1)/z evaluated at the argument of the exponential
    %     function in the integrand.
    
    % Initialize
    L20 = Mesh.EdgeData(i).L2;      % old entries of L2 inner product
    L2 = zeros(nDofs);              % prealocated memory for new entries
    Dir = vertcat(Mesh.ElemData(ind).Dir); % the propagation directions on adjacent elements
    realDir = all(imag(Dir(:))==0); % check for real propagation direction
    omega = Options.Omega;          % Helmholtz wave number
    Edge = Mesh.Coordinates(Mesh.Edges(i,:),:); % coordinates of endpoints
    BmA = Edge(2,:) - Edge(1,:);    % vector along edge
    DirBmA = Dir*BmA';              % scalar product of pw. propagation directions and edge
    numNewInd = sum([LocData(ind).nNewDofs]); % number of new indices
    numInd = sum([LocData(ind).nInd]); % number of old indices to copy
    doNewInd = numNewInd>0;         % calculate new terms
    doMixedInd = doNewInd && numInd>0; % calculate mixed terms
    
    % Compute element-dependent terms
    expX0 = zeros(nDofs,1);         % values of the basis functions at the first endpoint of the edge
    I0 = zeros(1,numInd);           % indices of relevant entries in L20
    I1 = zeros(1,numInd);           % corresponding indices in L2
    NewInd = zeros(1,numNewInd);    % indices of new elements in L2
    offsetDir = 0;
    offsetI0Val = 0;
    offsetI1Val = 0;
    offsetI_Ind = 0;
    offsetNewVal = 0;
    offsetNewInd = 0;
    for j=ind
      % Initialize
      nDir = Mesh.ElemData(j).nDofs;
      ind0 = offsetDir+1:offsetDir+nDir;
      
      % Compute indices
      indI = offsetI_Ind+1:offsetI_Ind+LocData(j).nInd;
      indN = offsetNewInd+1:offsetNewInd+LocData(j).nNewDofs;
      I0(indI) = offsetI0Val + LocData(j).I0;
      I1(indI) = offsetI1Val + LocData(j).I1;
      NewInd(indN) = offsetNewVal + LocData(j).NewInd;
      
      % Compute factor for integral
      if(doNewInd)
        x0 = Mesh.Coordinates(Mesh.Elements(j,1),:);
        dx0 = Edge(1,:) - x0;
        expX0(ind0) = exp(complex(0,1)*omega*Dir(ind0,:)*dx0');
      end
      
      % Update indices
      offsetDir = offsetDir + nDir;
      offsetI0Val = offsetI0Val + LocData(j).nDofs0;
      offsetI1Val = offsetI1Val + LocData(j).nDofs1;
      offsetNewVal = offsetNewVal + LocData(j).nDofs1;
      offsetI_Ind = offsetI_Ind + LocData(j).nInd;
      offsetNewInd = offsetNewInd + LocData(j).nNewDofs;
    end
    
    % Copy old entries into new matrix
    L2(I1,I1) = L20(I0,I0);
    
    % Compute factor for NewInd-NewInd terms
    if(doNewInd)
      Mloc = transpose(expX0(NewInd,ones(1,numNewInd))) ...
        .*conj(expX0(NewInd,ones(1,numNewInd)));
      Mloc = Mloc*Mesh.EdgeData(i).Length;

      % Compute exponentials for NewInd-NewInd terms
      dDirBmA = DirBmA(NewInd,ones(1,numNewInd)).' - conj(DirBmA(NewInd,ones(1,numNewInd)));
      ind0 = dDirBmA~=0;              % indices where the integrand is not constant
      ind0u = triu(ind0,1);           % upper triangular part - lower triangular part is complex conjugate
      expBmA = zeros(numNewInd);
      if(realDir)
        x = complex(0,1)*omega*dDirBmA(ind0u);
        expBmA(ind0u) = expm1(x)./x;    % integral for upper triangular part
        expBmA = expBmA + expBmA';      % fill in lower triangular part
      else
        x = complex(0,1)*omega*dDirBmA(ind0);
        expBmA(ind0) = expm1(x)./x;
      end
      expBmA(~ind0) = 1;              % and one for constant terms
      Mloc = Mloc.*expBmA;            % update local mass matrix

      % Copy NewInd-NewInd terms into matrix
      L2(NewInd,NewInd) = Mloc;
    end
    
    % Compute factor for mixed terms
    if(doMixedInd)
      Mloc = transpose(expX0(NewInd,ones(1,numInd))) ...
        .*conj(expX0(I1,ones(1,numNewInd)));
      Mloc = Mloc*Mesh.EdgeData(i).Length;

      % Compute exponentials for mixed terms
      dDirBmA = DirBmA(NewInd,ones(1,numInd)).' - conj(DirBmA(I1,ones(1,numNewInd)));
      ind0 = dDirBmA~=0;              % indices where the integrand is not constant
      expBmA = zeros(numInd,numNewInd);
      x = complex(0,1)*omega*dDirBmA(ind0);
      expBmA(ind0) = expm1(x)./x;     % integral for nonconstant terms
      expBmA(~ind0) = 1;              % one for constant terms
      Mloc = Mloc.*expBmA;            % update local mass matrix

      % Copy mixed terms into matrix
      L2(I1,NewInd) = Mloc;
      L2(NewInd,I1) = Mloc';          % by symmetry
    end
    
    % Update data structure
    Mesh.EdgeData(i).L2 = L2;
    
  end
  
%% Generate output arguments
  
  % Reset superfluous arguments
  blacklist = {'Elements','RemDir','Dir','nDofs','DofData'};
  for i=1:length(blacklist)
    Options.(blacklist{i}) = [];
  end
  
  % Add options to mesh data structure
  Mesh.Data.Options = Options;
  
  % Do postcomputations on mesh
  if(isa(Options.Post,'function_handle'))
    Mesh = Options.Post(Mesh);
  elseif(isa(Options.Post,'cell'))
    for j=1:numel(Options.Post)
      Mesh = Options.Post{j}(Mesh);
    end
  end
  
  % Return mesh as first argument
  varargout{1} = Mesh;
  
  % Return various data as second argument
  if(nargout>=2)
    
    % Copy Options data structure
    Data.Options = Options;
    
    % Construct global basis transformation
    nInd = sum([LocData.nInd]);
    Data.I0 = zeros(nInd,1);
    Data.I1 = zeros(nInd,1);
    i0 = 0;
    i1 = 0;
    j = 0;
    for i=1:nElements
      ind = j+1:j+LocData(i).nInd;
      Data.I0(ind) = i0 + LocData(i).I0;
      Data.I1(ind) = i1 + LocData(i).I1;
      j = j + LocData(i).nInd;
      i0 = i0 + LocData(i).nDofs0;
      i1 = i1 + LocData(i).nDofs1;
    end
    
    varargout{2} = Data;
  end
  
return

%% Subfunctions

% Compute angles from direction vectors
function Phi = getPhi(Dir)

  nrmDir = sqrt(sum(real(Dir).^2,2));
  Dir = real(Dir)./nrmDir(:,[1 1]);
  
  Phi = acos(Dir(:,1));
  ind = Dir(:,2)<0;
  Phi(ind) = 2*pi - Phi(ind);

return