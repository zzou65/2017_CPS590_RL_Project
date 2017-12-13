% driver routine for several diffusion convection discretizations
% see below
   clear Mesh
  % Initialize constants
  JIG=2;
  d1=1;                      %impact of supg modification
  d2=1;
  a=10^-10               %amount of diffusivity
  NREFS =6;
  % select test code
  d=getData(5);
     
  QuadRule = P7O6(); 
  
  Mesh.Coordinates =d.Coordinates; 
  Mesh.Elements = d.Elements;
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = d.boundtype;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  err=zeros(NREFS,4);
  h=zeros(NREFS,1);
  Dofs=zeros(NREFS,1);
  
  for i = 1:NREFS
   
   %refine Mesh   
   Mesh = refine_REG(Mesh);
   
   % Mesh preprocessing
    
        switch(JIG)
        case 1
          NewMesh = Mesh;      
        case 2
          Loc = get_BdEdges(Mesh);
          Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
          FixedPos = zeros(size(Mesh.Coordinates,1),1);
          FixedPos(Loc) = 1;
          NewMesh = jiggle(Mesh,FixedPos);   
        case 3
          Loc = get_BdEdges(Mesh);
          Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
          FixedPos = zeros(size(Mesh.Coordinates,1),1);
          FixedPos(Loc) = 1;
          NewMesh = smooth(Mesh,FixedPos);
        end
   h(i)=get_MeshWidth(NewMesh);
   
   %Laplace
   A = assemMat_LFE(NewMesh,@STIMA_Lapl_LFE);
   %Mass
   M = assemMat_LFE(NewMesh,@MASS_LFE);
  
   %Convection term (upwinding)
   %diagonal Mass Matrix
   MassZero=assemMat_Mass0fD(NewMesh); 
   % upwinding Quadrature Rule= discrete differential forms
   D_up=assemMat_LFE(NewMesh, @STIMA_ContrGrad_Up, d.V_Handle);
%               % or equivalently
%               TopGrad=assemMat_TopGrad(NewMesh);
%               ContrOne=assemMat_Contr1f(NewMesh,d.V_Handle);
%               D=ContrOne*TopGrad;
   D_up=MassZero*D_up;          
%   D_up=M*D_up;
      
   % standard Galerkin
   D_st=assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,d.V_Handle,P1O2());
   
   % edge-Element interpolation of velocity
   % D_int=assemMat_LFE(NewMesh,@STIMA_ConvInt_LFE,d.V_Handle,gauleg(0,1,3));
   
   % SUPG
   % D_supg1 = assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,d.V_Handle,P1O2());
   D_supg2 = assemMat_LFE(NewMesh,@STIMA_SUPG_LFE,P1O2(), d.V_Handle,a,d1,d2);
   D_supg  = D_st+D_supg2;
   
   % self-cooked
   %M_QFELFE=assemMat_MASS_QFELFE(NewMesh,P3O3());
   %ContrOneQFE=assemMat_Contr1f_QFE(NewMesh,d.V_Handle);  % contraction of one forms
   %D_c=M_QFELFE*ContrOneQFE*TopGrad;
   %Laplace + convection
  
   A_up =(a*A+D_up);
   A_st =(a*A+D_st);
   %A_int =(a*A+D_int);
   A_supg =(a*A+D_supg);
 %  A_c=(a*A+D_c);
   % source term
  
   L = assemLoad_LFE(NewMesh,P1O2(),d.SOL_Handle,a);
   L_U = assemLoad_LFE(NewMesh,P1O2(),d.U_EX_Handle,a);
   L_supg=assemLoad_LFE_SUPG(NewMesh,P1O2(),d.V_Handle,d.SOL_Handle,a,d1,d2);
   
   % Direchlet boundary
  
   [U_up,FreeDofs] = assemDir_LFE(NewMesh,-1,d.U_EX_Handle,a);
   U_st=U_up;
  % U_int=U_up;
   U_supg=U_up;
   U_c=U_up;
   U=U_up;
   
   L_up = L - A_up*U_up;
   L_st = L - A_st*U_st;
  % L_c = L - A_c*U_c;
   L_supg = (L+L_supg) - A_supg*U_supg;
   
   % solving system
  
   U_up(FreeDofs) = A_up(FreeDofs,FreeDofs)\L_up(FreeDofs);
   U_st(FreeDofs) = A_st(FreeDofs,FreeDofs)\L_st(FreeDofs);
%   U_c(FreeDofs) = A_c(FreeDofs,FreeDofs)\L_c(FreeDofs);
   U_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);
   err(i,1) = L2Err_LFE(NewMesh,U_up,P7O6(),d.U_EX_Handle,0,a);
   err(i,2) = L2Err_LFE(NewMesh,U_st,P7O6(),d.U_EX_Handle,0,a);
   err(i,3) = L2Err_LFE(NewMesh,U_supg,P7O6(),d.U_EX_Handle,0,a)
  % err(i,4) = L2Err_LFE(NewMesh,U_c,P7O6(),d.U_EX_Handle,0,a)
   Dofs(i)=size(NewMesh.Coordinates,1);
   
   plot_LFE(U_up,NewMesh); colorbar;
   plot_LFE(U_st,NewMesh);  colorbar;
  plot_LFE(U_supg,NewMesh); colorbar;

  % plot_LFE(U_c,NewMesh); colorbar;
  end; 
  U(FreeDofs)=M(FreeDofs,FreeDofs)\L_U(FreeDofs);
  plot_LFE(U,NewMesh);  
  grid on;
  plot_Mesh(NewMesh,'a')  
  fig = figure('Name','Discretization error');
  plot(h,err(:,1),'r--',h,err(:,2),'b-.',h,err(:,3),'g-',h,err(:,4),'c-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('L^2-error u (Up)','L^2-error u','L^2-error u (SUPG)','Location','NorthEast')
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,1)),1);
  add_Slope(gca,'SouthEast',p(1),'r--');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,2)),1);
  add_Slope(gca,'East',p(1),'b-.');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,3)),1);
  add_Slope(gca,'South',p(1),'g-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,4)),1);
  add_Slope(gca,'Southwest',p(1),'c-');
clear all;
