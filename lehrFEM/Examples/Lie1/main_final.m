%

  % Initialize constants
  
  NREFS =6;
  F_HANDLE = @f_;       % Right hand side source term
  GD_HANDLE = @g_;      % Dirichlet boundary data
  V_HANDLE = @velo;      % Velocity data 
  
  
%  Initialize mesh
  
%   Mesh.Coordinates = [-1 -1; ...
%                        1 -1; ...
%                        1  1; ... 
%                       -1  1];
%   Mesh.Elements = [1 2 3; ...
%                    1 3 4];


%   rand('state',sum(100*clock));
%   x = [rand(1,10),0,1,1,0];
%   y = [rand(1,10),0,0,1,1];
%   TRI = delaunay(x,y);
%   Mesh.Coordinates = [x',y'];
%   Mesh.Elements = TRI;
%   orient_Elems(Mesh);
  
  Mesh.Coordinates = [0 0; ...
                       1 0; ...
                       1  1; ... 
                      0  1];
  Mesh.Elements = [1 2 3; ...
                   1 3 4];
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);

  for i = 1:NREFS
    Mesh = refine_REG(Mesh);
  end

  % Laplace 
  
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  
  %Convection term (upwinding)

  % diagonal Mass Matrix

  MassZero=assemMat_MassZeroD(Mesh);
    
  % and generalized finite differences    
    
  D_fd=assemMat_LFE(Mesh, @LIEUP_FD, V_HANDLE);  
  

  % Laplace + convection
  
  d=10^(-3);
  c=1;
  
  A_fd =(d*A+c*MassZero*D_fd);
  
  % source term
  
  L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
  
  % Direchlet boundary
  
  [U_fd,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
  
  L_fd = L - A_fd*U_fd;
   
  % solving system
  
  U_fd(FreeDofs) = A_fd(FreeDofs,FreeDofs)\L_fd(FreeDofs);
  
  % plot solution

  plot_LFE(U_fd,Mesh);
  colorbar;
  % Clear memory
  
 clear all;
 