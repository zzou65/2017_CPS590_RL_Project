function varargout = plot_hp(U,Mesh,Elem2Dof,p)
% PLOT_HP Plot piecewise polynomial functions on a mesh.
%
%   PLOT_HP(U,MESH,ELEM2DOF,P) generates a plot of the piecewise polynomial
%   function of maximum polynomial degree P given by U on the struct MESH. 
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%
%   The struct ELEM2DOF specifies the element to dof mapping.
%
%   Example:
%
%   plot_HP(U,Mesh,Elem2Dof,3);
%
%   See also shap_hp.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  Offset = 0.05;                            % Axes offset parameter
  nCoordinates = size(Mesh.Coordinates,1);  % Number of vertices
  nElements = size(Mesh.Elements,1);        % Number of elements
  
  % Split the reference element according to maximum polynomial degree

  [RefCoord,RefElem] = split_RefElem(1/(2*p));
  nRefCoord = size(RefCoord,1);
  nRefElem = size(RefElem,1);
  
  % Compute values of shape functions on the reference element
  
  Shap = shap_hp(RefCoord,p);
    
  % Preallocate memory
  
  Coordinates = zeros(nRefCoord*nElements,2);
  Elements = zeros(nRefElem*nElements,3); 
  FVal = zeros(nRefCoord*nElements,1);
  
  % Build virtual plotting mesh COORDINATES and ELEMENTS 
 
  CoordPtr = 0;
  CoordLoc = 1:nRefCoord;
  ElemPtr = 0;
  ElemLoc = 1:nRefElem;
  for i = 1:nElements
        
    % Compute element mapping
      
    bK = Mesh.Coordinates(Mesh.Elements(i,1),:);
    BK = [Mesh.Coordinates(Mesh.Elements(i,2),:)-bK; ...
          Mesh.Coordinates(Mesh.Elements(i,3),:)-bK];
     
    % Insert elements and coordinates
     
    Coordinates(CoordPtr+CoordLoc,:) = RefCoord*BK+ones(nRefCoord,1)*bK;
    Elements(ElemPtr+ElemLoc,:) = RefElem+CoordPtr;
      
    % Compute function value on current element
      
    vidx = Mesh.Elements(i,:);
    for j1 = 1:3
      FVal(CoordPtr+CoordLoc) = FVal(CoordPtr+CoordLoc) ...
                              + U(vidx(j1))*Shap.vshap.shap{j1};
      dofs = Elem2Dof.EDofs{j1}.Dofs{i};
      ndofs = Elem2Dof.EDofs{j1}.nDofs(i);
      dir = Elem2Dof.EDofs{j1}.Dir(i);
      for j2 = 1:ndofs
        FVal(CoordPtr+CoordLoc) = FVal(CoordPtr+CoordLoc) ...
                      + U(dofs(j2))*dir^(j2+1)*Shap.eshap{j1}.shap{j2};
      end
    end
    dofs = Elem2Dof.CDofs.Dofs{i};
    ndofs = Elem2Dof.CDofs.nDofs(i);
    for j1 = 1:ndofs
      FVal(CoordPtr+CoordLoc) = FVal(CoordPtr+CoordLoc) ...
                              + U(dofs(j1))*Shap.cshap.shap{j1};      
    end
    
    % Update pointers
    
    CoordPtr = CoordPtr+nRefCoord;
    ElemPtr = ElemPtr+nRefElem;      
  end
  
  % Generate plot
  
  XMin = min(Mesh.Coordinates(:,1));
  XMax = max(Mesh.Coordinates(:,1));
  YMin = min(Mesh.Coordinates(:,2));
  YMax = max(Mesh.Coordinates(:,2));
  CMin = min(FVal);
  CMax = max(FVal);
  
  XLim = [XMin XMax]+Offset*(XMax-XMin)*[-1 1]; 
  YLim = [YMin YMax]+Offset*(YMax-YMin)*[-1 1];
  CLim = [CMin CMax]+Offset*(CMax-CMin)*[-1 1];
  
  fig = figure('Name','hp finite elements');
  hold on;
  patch('Faces',Elements, ...
        'Vertices',Coordinates, ...
        'CData',FVal, ...
        'FaceColor','interp', ...
        'EdgeColor','none');
  patch('Faces',Mesh.Elements, ...
        'Vertices',Mesh.Coordinates, ...
        'FaceColor','none', ...
        'EdgeColor','k');
  hold off;
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

function [RefCoord,RefElem] = split_RefElem(res)
% SPLIT_REFELEM Splits the reference element.
%
%   [REFCOORD,REFELEM] = SPLIT_REFELEM(RES) splits the reference triangle
%   into smaller elements according to the resolution RES.
%
%   Example:
%
%   [RefCoord,RefElem] = split_RefElem(0.1);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
    
  % Generate triangulation of reference triangle  
      
  [x,y] = meshgrid(0:res:1);  
  loc = (y < 1+res/2-x);
  x = x(loc);
  y = y(loc);
  RefCoord = [x(:) y(:)];
  RefElem = delaunay(x,y);
    
  % Adjust orientation to counter-clockwise
    
  for i = 1:size(RefElem,1)
    AB = RefCoord(RefElem(i,2),:)-RefCoord(RefElem(i,1),:);
    AC = RefCoord(RefElem(i,3),:)-RefCoord(RefElem(i,1),:);
    if(AB(1)*AC(2)-AB(2)*AC(1) < 0)
      tmp = RefElem(i,1);
      RefElem(i,1) = RefElem(i,2);
      RefElem(i,2) = tmp;
    end
  end 
  
return
