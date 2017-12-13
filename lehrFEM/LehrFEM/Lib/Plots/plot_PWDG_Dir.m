function varargout = plot_PWDG_Dir(Mesh,c,normalize)
%PLOT_PWDG_DIR plot of plane wave basis functions
%   
%   PLOT_PWDG_DIR(MESH) plots the mesh represented by the data structure
%   MESH and the propagation directions of the discontinuous plane wave
%   basis functions.  These can be added to the data structure using
%   SET_DATA_PWDG.
%
%   PLOT_PWDG_DIR(MESH,C) scales the plane wave propagation directions
%   by the corresponding entries of the vector C.  C can, for example, be a
%   finite element solution.
%
%   PLOT_PWDG_DIR(MESH,C,NORMALIZE) normalizes the lengths of the vectors
%   representing propagation directions on each element if NORMALIZE is
%   true (default) and doesn't if NORMALIZE is false.
%
%   H = PLOT_PWDG_DIR(...) returns a handle to the figure.
%
%   See also set_Data_PWDG, plot_PWDG.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  nElements = size(Mesh.Elements,1);
  nVert = size(Mesh.Elements,2);
  nDofs = [Mesh.ElemData.nDofs];
  nDofsTot = sum(nDofs);
  OFFSET = 0.05;
  EDGECOLOR = 'b';
  VECCOLOR = 'k';

  % Set default arguments
  if(nargin<2 || isempty(c))
    c = ones(nDofsTot,1);
  else
    c = abs(c);
  end
  if(nargin<3 || isempty(normalize))
    normalize = true;
  end
    
  % Compute axes limits
  X = Mesh.Coordinates(:,1);
  Y = Mesh.Coordinates(:,2);
  XMin = min(X);
  XMax = max(X);
  YMin = min(Y);
  YMax = max(Y);
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1];
  YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
  
  % Compute midpoints of elements
  MidPts0 = zeros(nElements,2);
  for j=1:nVert
    MidPts0 = MidPts0 + Mesh.Coordinates(Mesh.Elements(:,j),:);
  end
  MidPts0 = MidPts0/nVert;
  
  % Construct vector containing mindpoints corresponding to each degree of
  % freedom
  MidPts = zeros(nDofsTot,2);
  for j=1:nElements
    MidPts(Mesh.ElemData(j).Ind,:) = MidPts0(j(ones(nDofs(j),1)),:);
  end
  
  % Normalize coefficients on each element
  if(normalize)
    for j=1:nElements
      ind = Mesh.ElemData(j).Ind;
      c(ind) = c(ind)/max(c(ind));
    end
  end
  
  % Construct global vector containing propagation directions of all plane
  % wave basis functions
  Dir = vertcat(Mesh.ElemData.Dir);
  Dir = Dir.*c(:,[1 1]);
  
  % Plot mesh and propagation directions
  fig = figure('Name','PWDG Mesh plot');
  quiver(MidPts(:,1),MidPts(:,2),Dir(:,1),Dir(:,2),1,VECCOLOR);
  patch('Faces', Mesh.Elements, ...
        'Vertices', Mesh.Coordinates, ...
        'FaceColor', 'none', ...
        'EdgeColor', EDGECOLOR);
  axis('equal');
  box('on');
  set(gca,'XLim',XLim,'YLim',YLim,'ZLim',[0 1],'CLim',[0 1]);
  title('\bf Plane wave basis functions');
  
  % Assign output arguments    
  if(nargout > 0)
    varargout{1} = fig;
  end
  
return