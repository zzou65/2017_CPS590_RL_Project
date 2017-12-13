%

  % Initialize constants
  
  NREFS =8;
  V_HANDLE = @velo;     % Velocity data 
  W_EX_Handle=@(x,varargin) sin(pi.*x(:,1))+cos(pi.*x(:,2));
  W_GRAD_Handle=@(x,varargin) [pi.*cos(pi.*x(:,1)) -pi.*sin(pi.*x(:,2))];
  W_SOL_Handle=@(x,varargin) pi.^2.*sin(pi.*x(:,1))+pi.^2.*cos(pi.*x(:,2))+pi.*cos(pi.*x(:,1))-pi.*sin(pi.*x(:,2));
  V_Handle=@(x,varargin)ones(size(x,1),2);
  
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
  err=zeros(NREFS,1);
  err2=zeros(NREFS,1);
  h=zeros(NREFS,1);
  Dofs=zeros(NREFS,1);
  for i = 1:NREFS
   
   % refine Mesh   
   Mesh = refine_REG(Mesh);

   % Laplace
   A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  
   %Convection term (upwinding)

    % diagonal Mass Matrix
   MassZero=assemMat_MassZeroD(Mesh);
    
   % and generalized finite differences    
 %  D_fd=assemMat_LFE(Mesh, @LIEUP_FD, V_Handle);  
  
   TopGrad=assemMat_TopGrad(Mesh);
   ContrOne=assemMat_ContrOne(Mesh,V_Handle);
   
   D=ContrOne*TopGrad;
   
   % Laplace + convection
  
   A_fd =(A+MassZero*D);
  
   % source term
  
   L = assemLoad_LFE(Mesh,P7O6(),W_SOL_Handle);
  
   % Direchlet boundary
  
   [U_fd,FreeDofs] = assemDir_LFE(Mesh,-1,W_EX_Handle);
  
   L_fd = L - A_fd*U_fd;
   
   % solving system
  
   U_fd(FreeDofs) = A_fd(FreeDofs,FreeDofs)\L_fd(FreeDofs);
  
   err(i) = L2Err_LFE(Mesh,U_fd,P7O6(),W_EX_Handle)
   err2(i) = H1SErr_LFE(Mesh,U_fd,P7O6(),W_GRAD_Handle);
   h(i)=get_MeshWidth(Mesh);
   Dofs(i)=size(Mesh.Coordinates,1);

  end; 
  
  fig = figure('Name','Discretization error');
  plot(Dofs,err,'ro--',Dofs,err2,'bo--'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf Dofs}');
  ylabel('{\bf Error}');
  legend('L^2-error u','H^1S-error u','Location','NorthEast')
  p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(NREFS-3:NREFS)),1);
  add_Slope(gca,'SouthEast',p(1),'r-');
  p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS)),1);
  add_Slope(gca,'East',p(1),'b-');
clear all;