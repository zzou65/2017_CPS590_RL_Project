%

  % Initialize constants
  
  e=10^-5;                                  %diffusivity    
  NREFS =6;
  V_HANDLE = @(x,varargin) 1/sqrt(10)*[1,3];     % Velocity data 
  W_EX_Handle=@(x,varargin) sin(pi.*x(:,1))+cos(pi.*x(:,2));
  W_GRAD_Handle=@(x,varargin) [pi.*cos(pi.*x(:,1)) -pi.*sin(pi.*x(:,2))];
  W_SOL_Handle=@(x,varargin) e*(pi.^2.*sin(pi.*x(:,1))+pi.^2.*cos(pi.*x(:,2)))+pi.*cos(pi.*x(:,1))-pi.*sin(pi.*x(:,2));
  V_Handle=@(x,varargin)ones(size(x,1),2);
    
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
  err0=zeros(NREFS,1);
  err1=zeros(NREFS,1);
  err2=zeros(NREFS,1);
  err3=zeros(NREFS,1);
  h=zeros(NREFS,1);
  Dofs=zeros(NREFS,1);
  for i = 1:NREFS
   
   %refine Mesh   
   Mesh = refine_REG(Mesh);

   %Laplace
   A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  
   %Convection term (upwinding)

   %diagonal Mass Matrix
   MassZero=assemMat_Mass0fD(Mesh);
    
   % upwinding Quadrature Rule= discrete differential forms
   D_up=assemMat_LFE(Mesh, @STIMA_ContrGrad_Up, V_Handle);
   % or equivalently
   %TopGrad=assemMat_TopGrad(Mesh);
   %ContrOne=assemMat_Contr1f(Mesh,V_Handle);
   %D=ContrOne*TopGrad;
   
   
   % standard Galerkin
   D_st=assemMat_LFE(Mesh, @STIMA_ContrGrad_LFE, V_Handle);
   
   % Laplace + convection
  
   A_up =(e*A+MassZero*D_up);
   A_st =(e*A+D_st);
   % source term
  
   L = assemLoad_LFE(Mesh,P7O6(),W_SOL_Handle);
  
   % Direchlet boundary
  
   [U_up,FreeDofs] = assemDir_LFE(Mesh,-1,W_EX_Handle);
   U_st=U_up;
   
   L_up = L - A_up*U_up;
   L_st = L - A_st*U_st;
   
   % solving system
  
   U_up(FreeDofs) = A_up(FreeDofs,FreeDofs)\L_up(FreeDofs);
   U_st(FreeDofs) = A_st(FreeDofs,FreeDofs)\L_st(FreeDofs);
   err0(i) = L2Err_LFE(Mesh,U_up,P7O6(),W_EX_Handle)
   err1(i) = L2Err_LFE(Mesh,U_st,P7O6(),W_EX_Handle)
   err2(i) = H1SErr_LFE(Mesh,U_up,P7O6(),W_GRAD_Handle);
   err3(i) = H1SErr_LFE(Mesh,U_st,P7O6(),W_GRAD_Handle);
   h(i)=get_MeshWidth(Mesh);
   Dofs(i)=size(Mesh.Coordinates,1);
   
   plot_LFE(U_up,Mesh);
   plot_LFE(U_st,Mesh);
  end; 
  
  fig = figure('Name','Discretization error');
  plot(Dofs,err0,'ro--',Dofs,err1,'yo--',Dofs,err2,'go--',Dofs,err3,'bo--'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf Dofs}');
  ylabel('{\bf Error}');
  legend('L^2-error u (up)','L^2-error u','H^1S-error u (up)','H^1S-error u','Location','NorthEast')
  p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err0(NREFS-3:NREFS)),1);
  add_Slope(gca,'SouthEast',p(1),'r-');
  p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS)),1);
  add_Slope(gca,'East',p(1),'y-');
clear all;