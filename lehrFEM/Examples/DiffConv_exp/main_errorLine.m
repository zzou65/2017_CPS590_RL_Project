% driver routine for several diffusion convection discretizations
% see below
   clear Mesh
  % Initialize constants
  JIG=2;
  d1=0.5;                   %impact of supg modification
  d2=0.17;
  a=10^-10;               %amount of diffusivity
  NREFS =7;
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
  err1=zeros(NREFS,4);
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
   D=assemMat_LFE(NewMesh, @STIMA_ContrGrad_Up, d.V_Handle);
%               % or equivalently
%   TopGrad=assemMat_TopGrad(NewMesh);
%               ContrOne=assemMat_ContrOne(NewMesh,d.V_Handle);
%               D=ContrOne*TopGrad;
   D_up=MassZero*D;          
   D_up1=M*D;
      
   % standard Galerkin
   D_st=assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,d.V_Handle,P1O2());
   
   % edge-Element interpolation of velocity
   % D_int=assemMat_LFE(NewMesh,@STIMA_ConvInt_LFE,d.V_Handle,gauleg(0,1,3));
   
   % SUPG
   % D_supg1 = assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,d.V_Handle,P1O2());
   D_supg2 = assemMat_LFE(NewMesh,@STIMA_SUPG_LFE,P1O2(), d.V_Handle,a,d1,d2);
   D_supg  = D_st+D_supg2;
   
   % Upwinfd using midpoint quadrature
   UP=assemMat_UpLFE(NewMesh,d.V_Handle);
   
   
   % self-cooked
   % M_QFELFE=assemMat_MASS_QFELFE(NewMesh,P3O3());
   % ContrOneQFE=assemMat_Contr1f_QFE(NewMesh,d.V_Handle);  % contraction of one forms
   % D_c=M_QFELFE*ContrOneQFE*TopGrad;
   % Laplace + convection
  
   A_up =(a*A+D_up);
   A_up1 =(a*A+D_up1);
   A_up2 =(a*A+UP);
   A_st =(a*A+D_st);
   A_supg =(a*A+D_supg);
   % source term
  
   L = assemLoad_LFE(NewMesh,P7O6(),d.SOL_Handle,a);
   L_U = assemLoad_LFE(NewMesh,P7O6(),d.U_EX_Handle,a);
   L_supg=assemLoad_LFE_SUPG(NewMesh,P7O6(),d.V_Handle,d.SOL_Handle,a,d1,d2);
   
   % Direchlet boundary
  
   [U_up,FreeDofs] = assemDir_LFE(NewMesh,-1,d.U_EX_Handle,a);
   U_st=U_up;
   U_up1=U_up;
   U_up2=U_up;
   U_supg=U_up;
   
   L_up = L - A_up*U_up;
   L_st = L - A_st*U_st;
   L_up1 = L - A_up1*U_up1;
   L_up2 = L - A_up2*U_up2;
   L_supg = (L+L_supg) - A_supg*U_supg;
   
   % solving system
  
   U_up(FreeDofs) = A_up(FreeDofs,FreeDofs)\L_up(FreeDofs);
   U_st(FreeDofs) = A_st(FreeDofs,FreeDofs)\L_st(FreeDofs);
   U_up1(FreeDofs) = A_up1(FreeDofs,FreeDofs)\L_up1(FreeDofs);
   U_up2(FreeDofs) = A_up2(FreeDofs,FreeDofs)\L_up2(FreeDofs);
   U_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);
   err(i,1) = L2Err_LFE(NewMesh,U_up,P7O6(),d.U_EX_Handle,0,a);
   err(i,2) = L2Err_LFE(NewMesh,U_st,P7O6(),d.U_EX_Handle,0,a);
   err(i,3) = L2Err_LFE(NewMesh,U_supg,P7O6(),d.U_EX_Handle,0,a);
   err(i,4) = L2Err_LFE(NewMesh,U_up1,P7O6(),d.U_EX_Handle,0,a);
   err(i,5) = L2Err_LFE(NewMesh,U_up2,P7O6(),d.U_EX_Handle,0,a)
   err1(i,1) = L2Err_LFEmod(NewMesh,U_up,P7O6(),d.U_EX_Handle,0,a);
   err1(i,2) = L2Err_LFEmod(NewMesh,U_st,P7O6(),d.U_EX_Handle,0,a);
   err1(i,3) = L2Err_LFEmod(NewMesh,U_supg,P7O6(),d.U_EX_Handle,0,a);
   err1(i,4) = L2Err_LFEmod(NewMesh,U_up1,P7O6(),d.U_EX_Handle,0,a);
   err1(i,5) = L2Err_LFEmod(NewMesh,U_up2,P7O6(),d.U_EX_Handle,0,a);
   Dofs(i)=size(NewMesh.Coordinates,1);
   
   plot_LFE(U_up,NewMesh); colorbar;
   plot_LFE(U_up2,NewMesh);  colorbar;
   plot_LFE(U_supg,NewMesh); colorbar;

   %plot_LFE(U_c,NewMesh); colorbar;
   %plot_Mesh(NewMesh,'a')
  end;
  x0=[0,0];x1=[1,1];
  [e1, u1]=plotLine_LFE(U_up,NewMesh,x0, x1); colorbar;
  [e2, u2]=plotLine_LFE(U_supg,NewMesh,x0, x1); colorbar;
  [e3, u3]=plotLine_LFE(U_up1,NewMesh,x0, x1); colorbar;
  [e4, u4]=plotLine_LFE(U_up2,NewMesh,x0, x1); colorbar;
  e5=e4; u5=d.U_EX_Handle(ones(size(e4,2),1)*x0+e4'*(x1-x0)/norm(x1-x0),0,a)';
  
  Xmin=min([min(e1),min(e2),min(e3),min(e4),min(e5)]); Xmax=max([max(e1),max(e2),max(e3),max(e4),max(e5)]);
  Ymin=min([min(u1),min(u2),min(u3),min(u4),min(u5)]); Ymax=max([max(u1),max(u2),max(u3),max(u4),max(u5)]);
  figure;
  plot(e1,u1,'r-',e2,u2,'b-',e3,u3,'r--',e4,u4,'g-',e5,u5,'b--');
  set(gca,'Xlim',[Xmin-0.05,Xmax+0.05],'Ylim',[Ymin-0.05,Ymax+0.05])
  legend('upwind (lumped)', 'SUPG','upwind (full)','upwind (midpoints)','Location','NorthWest');
  xlabel('{\bf profile line}');ylabel('{\bf profile }');
 % print -depsc ~/smoothing_profile.eps
 % plot_LFE(U_up,NewMesh); colorbar;
 % print -depsc ~/smoothing_upwind.eps
 % plot_LFE(U_supg,NewMesh); colorbar;
 % print -depsc ~/smoothing_supg.eps
  h(end);

  % plot_LFE(M\L_U,NewMesh); colorbar;  
  fig = figure('Name','Discretization error');
  plot(h,err(:,1),'r--',h,err(:,2),'b-.',h,err(:,3),'g-',h,err(:,4),'c-',h,err(:,5),'y-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('L^2-error u (Up-lumped)','L^2-error u','L^2-error u (SUPG)','L^2-error u (Up-full)','L^2-error u (midpoints)','Location','NorthWest')
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,1)),1);
  add_Slope(gca,'SouthEast',p(1),'r--');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,2)),1);
  add_Slope(gca,'East',p(1),'b-.');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,3)),1);
  add_Slope(gca,'NorthEast',p(1),'g-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,4)),1);
  add_Slope(gca,'Southwest',p(1),'c-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,5)),1);
  add_Slope(gca,'Southwest',p(1),'y-');
%  print -depsc converge_compl_-10.eps
%   print -dfig converge_compl_-10.fig
%    
%   fig = figure('Name','Discretization error');
%   plot(h,err1(:,1),'r--',h,err1(:,2),'b-.',h,err1(:,3),'g-',h,err1(:,4),'c-'); grid('on');
%   set(gca,'XScale','log','YScale','log');
%   xlabel('{\bf h}');
%   ylabel('{\bf Error}');
%   legend('L^2-error u (Up)','L^2-error u','L^2-error u (SUPG)','Location','NorthWest')
%   p = polyfit(log(h(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS,1)),1);
%   add_Slope(gca,'SouthEast',p(1),'r--');
%   p = polyfit(log(h(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS,2)),1);
%   add_Slope(gca,'East',p(1),'b-.');
%   p = polyfit(log(h(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS,3)),1);
%   add_Slope(gca,'NorthEast',p(1),'g-');
%   p = polyfit(log(h(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS,4)),1);
%   add_Slope(gca,'Southwest',p(1),'c-');
%   print -depsc converge_part_-10.eps
clear all;