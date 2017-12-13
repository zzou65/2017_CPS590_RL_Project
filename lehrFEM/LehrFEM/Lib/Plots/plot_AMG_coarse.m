function varargout = plot_AMG_coarse(Mesh,AMGData,FDofs)
%PLOT_AMG_COARSE plot coarse grid points for amg
%
%   PLOT_AMG_COARSE(MESH,AMGDATA) plots the mesh MESH and highlights the
%   coarse grid points corresponding to the AMG data structure AMGDATA.
%   AMGDATA can be created with AMGSETUP or SOLVE_AMG.
%
%   PLOT_AMG_COARSE(MESH,AMGDATA,FDOFS) only plots points in
%   MESH.Coordinates(FDOFS,:).  This form is required if AMGDATA does not
%   handle all degrees of freedom.
%
%   FIG = PLOT_AMG_COARSE(...) returns the handle to the figure.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the
%                mesh. 
%
%   Example:
%
%       plot_AMG_coarse(Mesh,AMGSetup(A));
%
%   See also plot_Mesh, amg_solve.


%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % If FDofs argument is given, use only corresponding poins
  
  if(nargin >= 3 && ~isempty(FDofs))
    Coordinates = Mesh.Coordinates(FDofs,:);
  else
    Coordinates = Mesh.Coordinates;
  end

  %% Determine level of grid points

  % maximal level

  LVL = length(AMGData);
  
  % recursively calculate level of all points
  
  level = ones(size(AMGData{LVL}.A,1),1);
  
  for j=2:LVL
    levelc = level;
    level = j*ones(size(AMGData{LVL-j+1}.A,1),1);
    level(AMGData{LVL-j+1}.c) = levelc;
  end
  
  %% Plot coarse grid points
  
  % Determine size of circles
  
%   s = 5*(2.^((0:LVL-1)/(LVL-1)));
  s = 10*(LVL+1-(0:LVL))/(LVL+1);

  % Plot Mesh
  
  fig = figure('Name','AMG Coarse Grid Points');
  plot_Mesh(Mesh,'af');
  
  % Initialize plot for levels
  
  hold all;
  lgd = cell(0);
  h = [];
  
  % plot levels
  
  for j=LVL-1:-1:1
    pts = (level==j);
    h1 = plot(Coordinates(pts,1),Coordinates(pts,2),'ro',...
      'MarkerSize',s(j),'MarkerFaceColor','r');
    lgd = {lgd{:},sprintf('level %.0f',j)};
    h = [h,h1];
  end
  
  if(length(lgd)>1)
    legend(h,lgd{:},'Location','NorthEastOutside');
  end
  
  title('\bf Coarse Grid Points');
  
  % Assign output arguments    
  
  switch(nargout)
    case 1
      varargout{1} = fig;
    case 2
      varargout{1} = fig;
      varargout{2} = h;
  end
return  