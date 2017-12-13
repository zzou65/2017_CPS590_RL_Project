function varargout = plot_PDG(u,Mesh,p,Shap)
% PLOT_PDG Plot discontinuous, piecewise polynomial functions on a mesh.
%
%   PLOT_PDG(U,MESH,P,SHAP) generates a plot of the discontinuous,
%   piecewise polynomial function of polynomial degree P given by U and the
%   reference element shape functions SHAP on the struct MESH. 
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%
%   SHAP is a function handle to the reference element shape functions.
%
%   Example:
%
%   plot_PDG(U,Mesh,2,shap_QFE);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  Offset = 0.05;                            % Offset parameter
  nCoordinates = size(Mesh.Coordinates,1);  % Number of vertex coordinates
  nElements = size(Mesh.Elements,1);        % Number of elements

  % Split the reference element

  [RefCoord,RefElem] = split_RefElem(p);
  nRefCoord = size(RefCoord,1);
  nRefElem = size(RefElem,1);
  
  % Evaluate shape functions
  
  N = Shap(RefCoord);
  nDofs = size(N,2);
  
  % Allocate memory
  
  FVal = zeros(nRefCoord*nElements,1);
  Coordinates = zeros(nRefCoord*nElements,2);
  Elements = zeros(nRefElem*nElements,3);
  
  CoordLoc = 1:nRefCoord;
  CoordPtr = 0;
  ElemLoc = 1:nRefElem;
  FLoc = 1:nDofs;
  for i = 1:nElements
    
    % Compute element mapping
    
    bK = Mesh.Coordinates(Mesh.Elements(i,1),:);
    BK = [Mesh.Coordinates(Mesh.Elements(i,2),:) - bK; ...
          Mesh.Coordinates(Mesh.Elements(i,3),:) - bK];
    
    % Assign auxiliary elements and coordinates 
    
    Coordinates(CoordLoc,:) = RefCoord*BK+ones(nRefCoord,1)*bK;
    Elements(ElemLoc,:) = RefElem + CoordPtr;
    
    % Assign function values
    
    FVal(CoordLoc) = N*u(FLoc);
    
    % Update locators
    
    CoordLoc = CoordLoc + nRefCoord;
    CoordPtr = CoordPtr + nRefCoord;
    ElemLoc = ElemLoc + nRefElem;
    FLoc = FLoc + nDofs;
    
  end
  
  % Generate plot
  
  XMin = min(Mesh.Coordinates(:,1));
  XMax = max(Mesh.Coordinates(:,1));
  YMin = min(Mesh.Coordinates(:,2));
  YMax = max(Mesh.Coordinates(:,2));
  CMin = min(FVal);
  CMax = max(FVal);
  
  XLim = [XMin XMax] + Offset*(XMax-XMin)*[-1 1]; 
  YLim = [YMin YMax] + Offset*(YMax-YMin)*[-1 1];
  CLim = [CMin CMax] + Offset*(CMax-CMin)*[-1 1];
  
  fig = figure('Name','hp finite elements');
  patch('Faces',Elements, ...
        'Vertices',[Coordinates FVal], ...
        'CData',FVal, ...
        'FaceColor','interp', ...
        'EdgeColor','none');
  colormap('jet');
  set(gca,'XLim',XLim, ...
          'YLim',YLim, ...
          'CLim',CLim, ...
          'DataAspectRatio',[1 1 1]);
  set(gcf,'Renderer','zbuffer');
 
  % Assign output argument
  
  if(nargout > 0)
    varargout{1} = fig;
  end
  
return

%%% Split reference element %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RefCoord,RefElem] = split_RefElem(nSteps)
% SPLIT_REFELEM Splits the reference element.
%
%   [REFCOORD,REFELEM] = SPLIT_REFELEM(NSTEPS) splits the reference
%   triangle into smaller elements according to the resolution NSTEPS.
%
%   Example:
%
%   [RefCoord,RefElem] = split_RefElem(3);
%
%   See also add_Edges, refine_REG.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
    
  % Initialize mesh
  
  Mesh.Coordinates = [0 0; 1 0; 0 1];
  Mesh.Elements = [1 2 3];
  Mesh = add_Edges(Mesh);
  Mesh.BdFlags = -ones(3,1);
  
  % Do regular mesh refinements
  
  for i = 1:nSteps
    Mesh = refine_REG(Mesh);
  end
  
  % Extract auxiliary elements and coordinates
  
  RefCoord = Mesh.Coordinates;
  RefElem = Mesh.Elements;
  
return