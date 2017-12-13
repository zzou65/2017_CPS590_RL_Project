% Obstackle 2D
%   Run script for obstacle problem on the unit square with zero boundary
%   conditions. The obstacle is given by 
%   Obstacle = (1-x(:,1).^2).*(1-x(:,2).^2)
%   Plots of the solution and the obstacle is saved as solution.eps and
%   obstacle.eps

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland


% Initialize constants


clc
warning off
TOL = 1e-10;
MAXIT = inf;
NREFS = 5;
isw = 2; % 1=PSOR, 2 = CG, 3 = quadprog
Omega = 2/(1+sqrt(2));
GD_HANDLE = @(x,varargin)0;
F_HANDLE = @(x,varargin)-50;
Obst = @(x,varargin)5*(1-x(:,1).^2).*(1-x(:,2).^2);
%Obst = @(x,varargin)1/8+1/2*sqrt(x(:,1).^2+x(:,2).^2)-3/4*(x(:,1).^2+x(:,2).^2);

% Initialize mesh

Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

for i = 1:NREFS,
    Mesh = refine_REG(Mesh);
end
Mesh = add_Edge2Elem(Mesh);

% Assemble stiffness matrix and load vector

A = assemMat_LFE(Mesh, @STIMA_Lapl_LFE);
L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);

% Assign obstacle

nCoord = size(Mesh.Coordinates,1);
psi = zeros(nCoord,1);
for i = 1:nCoord,
    coord = Mesh.Coordinates(i,:);
    psi(i) = Obst(coord);
end

% Incorporate boundary data (Dirichlet)

[U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
L = L - A*U;

% Shift problem

x = U(FreeDofs) - psi(FreeDofs);
b = L(FreeDofs) - A(FreeDofs,FreeDofs)*psi(FreeDofs);

% Calculate solution

if isw==1,
    [x conv] = psor(A(FreeDofs,FreeDofs),b,x,Omega,TOL,MAXIT);
elseif isw==2,
    [x conv] = cg(A(FreeDofs,FreeDofs),b,x,TOL,MAXIT);
elseif isw==3
    x = quadprog(A(FreeDofs,FreeDofs),-b,[],[],[],[],zeros(size(FreeDofs,2),1));
    conv = 1;
end
warning on

if ~conv,
    error('The system does not converge within the specified number of iterations')
end

U(FreeDofs) = x + psi(FreeDofs);

for i = 1:nCoord,
    if U(i)<psi(i),
        warning('Solution is not respecting the obstacle')
    end
end

% plot solution

  % Initialize constants
  
  OFFSET = 0.05;
  
  % Compute axes limits
  
  XMin = min(Mesh.Coordinates(:,1));
  XMax = max(Mesh.Coordinates(:,1));
  YMin = min(Mesh.Coordinates(:,2));
  YMax = max(Mesh.Coordinates(:,2));
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1]; 
  YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
  
  % Generate figure
    
  
    % Compute color axes limits 
      
    CMin = min(U);
    CMax = max(U);
    if(CMin < CMax)          % or error will occur in set function
      CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    else
      CLim = [1-OFFSET 1+OFFSET]*CMin;   
    end
         
    % Plot real finite element solution  
    
    fig = figure('Name', 'Solution');
    patch('faces', Mesh.Elements, ...
          'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) U], ...
          'CData', U, ...
          'facecolor', [.7 .4 .4], ...
          'edgecolor', 'none', ...
          'Facelighting', 'phong', ...
          'Ambientstrength', .3, ...
          'DiffuseStrength', .9, ...
          'Clipping', 'off', ...
          'BackFaceLighting', 'lit', ...
          'SpecularStrength', .9, ...
          'SpecularColorReflectance', 1, ...
          'SpecularExponent', 7);
    l1 = light('Position', [1 -.1 1.3], ...
               'Style', 'infinite', ...
               'Color', [1 1 0]);
    l2 = light('Position', [-1 -1 2], ...
               'Style', 'infinite', ...
               'Color', [.5 .4 0]);
    set(gca,'XLim',XLim, ...
            'YLim',YLim, ...
            'CLim',CLim, ...
            'DataAspectRatio',[1 1 2.5], ...
            'CameraPosition', [-9.1314 -11.9003 22.1506], ...
            'CameraUpVector', [0 0 1], ...
            'Visible', 'on', ...
            'LineWidth', 1.5);
    grid on
    
print -depsc solution.eps

% Plot Obstacle
     
  
    % Compute color axes limits 
      
    CMin = min(psi);
    CMax = max(psi);
    if(CMin < CMax)          % or error will occur in set function
      CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    else
      CLim = [1-OFFSET 1+OFFSET]*CMin;   
    end
         
    % Plot real finite element solution  
      
    fig = figure('Name','Obstacle');
    patch('faces', Mesh.Elements, ...
          'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) psi], ...
          'CData', psi, ...
          'facecolor', [.9 .2 .2], ...
          'edgecolor', 'none', ...
          'Facelighting', 'phong', ...
          'Ambientstrength', .3, ...
          'DiffuseStrength', .9, ...
          'Clipping', 'off', ...
          'BackFaceLighting', 'lit', ...
          'SpecularStrength', .9, ...
          'SpecularColorReflectance', 1, ...
          'SpecularExponent', 7);
    set(gca,'XLim',XLim, ...
            'YLim',YLim, ...
            'CLim',CLim, ...
            'DataAspectRatio',[1 1 2.5], ...
            'CameraPosition', [-9.1314 -11.9003 22.1506], ...
            'CameraUpVector', [0 0 1], ...
            'Visible', 'on', ...
            'LineWidth', 1.5);
    l1 = light('Position', [-1 -.5 1.3], ...
           'Style', 'infinite', ...
           'Color', [.8 .8 0]);
    l2 = light('Position', [1 1 2], ...
           'Style', 'infinite', ...
           'Color', [.5 .4 0]);
    grid on
print -depsc obstacle.eps