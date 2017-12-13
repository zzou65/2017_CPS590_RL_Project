function varargout = contour_PWDG(u,Mesh,doplots,varargin)
%CONTOUR_PWDG contour plot for plane wave DG
%
%   CONTOUR_PWDG(U,MESH) generates a contour plot of the real part of the
%   function represented by the vector U in the discontinuous plane wave
%   basis on the mesh MESH.  The plane wave information can be added to
%   MESH using SET_DATA_PWDG.
%
%   CONTOUR_PWDG(U,MESH,DOPLOTS) where DOPLOTS is a logical vector of
%   length three generates each of the following contour plots if the
%   corresponding entry of DOPLOTS is true (or 1):
%     1) real part
%     2) imaginary part
%     3) absolute value
%   If DOPLOTS is emtpy, only the real part is plotted.
%
%   CONTOUR_PWDG(U,MESH,DOPLOTS,PARAMS) passes the variable length
%   parameter list PARAMS to the plotting routine CONTOUR.
%
%   See also contour, eval_PWDG, plot_PWDG, plot_PWDG_Dir, set_Data_PWDG.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Set default arguments
  if(nargin<3 || isempty(doplots))
    doplots = [true false false];
  elseif(isa(doplots,'numeric'))
    doplots = logical(doplots);
  end

  % Initialize constants
  omega = Mesh.Data.Options.Omega;
  BDCOLOR = 'k';
  OFFSET = 0.05;

  % Initialize evaluation grid
  xmin = min(Mesh.Coordinates(:,1));
  xmax = max(Mesh.Coordinates(:,1));
  xlen = xmax - xmin;
  XLim = [xmin xmax] + OFFSET*xlen*[-1 1];

  ymin = min(Mesh.Coordinates(:,2));
  ymax = max(Mesh.Coordinates(:,2));
  ylen = ymax - ymin;
  YLim = [ymin ymax] + OFFSET*ylen*[-1 1];
  
  diam = sqrt(xlen^2 + ylen^2);
  h = max(2*pi/(50*omega),diam/1000);
  hx = xlen/ceil(xlen/h);
  hy = ylen/ceil(ylen/h);
  
  % Generate evaluation grid
  x = xmin:hx:xmax;
  y = ymin:hy:ymax;
  [X,Y] = meshgrid(x,y);
  Z = zeros(size(X));
  
  % Construct list of evaluation points
  XY = zeros(numel(X),2);
  XY(:,1) = X(:);
  XY(:,2) = Y(:);
  
  % Evaluate function u at gridpoints
  Z(:) = eval_PWDG(Mesh,omega,u,XY);
  
  % Find boundary edges for linear and parametrizend boundaries
  Loc = get_BdEdges(Mesh);
  areParam = isfield(Mesh,'EdgeData') && isfield (Mesh.EdgeData,'Geom');
  if(areParam)
    isParam = false(size(Loc));
    for j=1:length(Loc)
      isParam(j) = Mesh.EdgeData(Loc(j)).Geom.isParam;
    end
%     isParam = vertcat(Mesh.EdgeData(Loc').Geom.isParam);
    LocParam = Loc(isParam);
    Loc = Loc(~isParam);
  else
    LocParam = zeros(0,1);
  end  
  
  % Compute boundary edges for piecewise linear boundaries
%   Loc = get_BdEdges(Mesh);
  BdEdges_x = zeros(2,size(Loc,1));
  BdEdges_y = zeros(2,size(Loc,1));
  BdEdges_x(1,:) = Mesh.Coordinates(Mesh.Edges(Loc,1),1)';
  BdEdges_x(2,:) = Mesh.Coordinates(Mesh.Edges(Loc,2),1)';
  BdEdges_y(1,:) = Mesh.Coordinates(Mesh.Edges(Loc,1),2)';
  BdEdges_y(2,:) = Mesh.Coordinates(Mesh.Edges(Loc,2),2)';
  
  % Compute linear boundary edges for parametrized boundaries 
  BdParamEdges_x = zeros(2,size(LocParam,1));
  BdParamEdges_y = zeros(2,size(LocParam,1));
  BdParamEdges_x(1,:) = Mesh.Coordinates(Mesh.Edges(LocParam,1),1)';
  BdParamEdges_x(2,:) = Mesh.Coordinates(Mesh.Edges(LocParam,2),1)';
  BdParamEdges_y(1,:) = Mesh.Coordinates(Mesh.Edges(LocParam,1),2)';
  BdParamEdges_y(2,:) = Mesh.Coordinates(Mesh.Edges(LocParam,2),2)'; 
  
  % Generate plots
  fig = zeros(1,3);
  
  if(doplots(1)) % real part
    fig(1) = figure('Name','Plane wave DG finite elements');
    plot(BdEdges_x,BdEdges_y,BDCOLOR);
    hold on;
    plot(BdParamEdges_x,BdParamEdges_y,'Color',[1 1 1],'LineWidth',1); % erase parametrized boundary edges
    t = linspace(0,1,20); % parametrized boundary edges
    for j=LocParam'
      x = Mesh.EdgeData(j).Geom.Gamma(t);
      plot(x(:,1),x(:,2),BDCOLOR);
    end
    contour(X,Y,real(Z),varargin{:});
    hold off;
    axis equal;
    set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);
    title('\bf Real Part');
  end
  
  if(doplots(2)) % imaginary part
    fig(2) = figure('Name','Plane wave DG finite elements');
    plot(BdEdges_x,BdEdges_y,BDCOLOR);
    hold on;
    plot(BdParamEdges_x,BdParamEdges_y,'Color',[1 1 1],'LineWidth',1); % erase parametrized boundary edges
    t = linspace(0,1,20); % parametrized boundary edges
    for j=LocParam'
      x = Mesh.EdgeData(j).Geom.Gamma(t);
      plot(x(:,1),x(:,2),BDCOLOR);
    end
    contour(X,Y,imag(Z),varargin{:});
    hold off;
    axis equal;
    set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);
    title('\bf Imaginary Part');
  end
  
  if(doplots(3)) % real part
    fig(3) = figure('Name','Plane wave DG finite elements');
    plot(BdEdges_x,BdEdges_y,BDCOLOR);
    hold on;
    plot(BdParamEdges_x,BdParamEdges_y,'Color',[1 1 1],'LineWidth',1); % erase parametrized boundary edges
    t = linspace(0,1,20); % parametrized boundary edges
    for j=LocParam'
      x = Mesh.EdgeData(j).Geom.Gamma(t);
      plot(x(:,1),x(:,2),BDCOLOR);
    end
    contour(X,Y,abs(Z),varargin{:});
    hold off;
    axis equal;
    set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);
    title('\bf Absolute Value');
  end

  if(nargout > 0)
    varargout = num2cell(fig(doplots));
  end
  
return  