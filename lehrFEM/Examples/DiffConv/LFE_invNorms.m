% driver routine for several diffusion convection discretizations
% see below
   clear Mesh
  % Initialize constants
  JIG=2;
  d1=1;                      %impact of supg modification
  d2=1;
  aa=10.^-[0:12];               %amount of diffusivity
  na=size(aa,2);
  NREFS =3;
  % select test code
  d=getData(6);
     
  QuadRule = P7O6(); 
  
  Mesh.Coordinates =d.Coordinates; 
  Mesh.Elements = d.Elements;
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = d.boundtype;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  err=zeros(na,4);
  err1=zeros(na,4);
  h=zeros(na,1);
  Dofs=zeros(na,1);
  
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
  end
  
  for ia=1:na
   
   h(ia)=aa(ia);
   a=aa(ia);
   
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
%   TopGrad=assemMat_TopGrad(NewMesh);
%               ContrOne=assemMat_ContrOne(NewMesh,d.V_Handle);
%               D=ContrOne*TopGrad;
   D_up=MassZero*D_up;          
   %D_up=M*D_up;
      
   % standard Galerkin
   D_st=assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,d.V_Handle,P1O2());
   
   % edge-Element interpolation of velocity
   % D_int=assemMat_LFE(NewMesh,@STIMA_ConvInt_LFE,d.V_Handle,gauleg(0,1,3));
   
   % SUPG
   % D_supg1 = assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,d.V_Handle,P1O2());
   D_supg2 = assemMat_LFE(NewMesh,@STIMA_SUPG_LFE,P1O2(), d.V_Handle,a,d1,d2);
   D_supg  = D_st+D_supg2;
   
   A_up =(a*A+D_up);
   A_st =(a*A+D_st);
   A_supg =(a*A+D_supg);
   % source term
  
  % L = assemLoad_LFE(NewMesh,P1O2(),d.SOL_Handle,a);
  % L_U = assemLoad_LFE(NewMesh,P1O2(),d.U_EX_Handle,a);
  % L_supg=assemLoad_LFE_SUPG(NewMesh,P1O2(),d.V_Handle,d.SOL_Handle,a,d1,d2);
   
   % Direchlet boundary
  
   [U_up,FreeDofs] = assemDir_LFE(NewMesh,[-1 -2],d.U_EX_Handle,a);
    U_st=U_up;
    U_supg=U_up;
   
  % L_up = L - A_up*U_up;
  % L_st = L - A_st*U_st;
  % L_supg = (L+L_supg) - A_supg*U_supg;
   
   % solving system
  I_up=inv(full(A_up(FreeDofs,FreeDofs)));
  I_st=inv(full(A_st(FreeDofs,FreeDofs)));
  I_supg=inv(full(A_supg(FreeDofs,FreeDofs)));
  % U_up(FreeDofs) = A_up(FreeDofs,FreeDofs)\L_up(FreeDofs);
  % U_st(FreeDofs) = A_st(FreeDofs,FreeDofs)\L_st(FreeDofs);
  % U_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);
   err(ia,1) = abs(min(min(min(I_up,0))));
   err(ia,2) = abs(min(min(min(I_st,0))));
   err(ia,3) = abs(min(min(min(I_supg,0))))
   err1(ia,1) = norm(I_up,inf);
   err1(ia,2) = norm(I_st,inf);
   err1(ia,3) = norm(I_supg,inf);
  end; 
  % plot_LFE(M\L_U,NewMesh); colorbar;  
  fig = figure('Name','Discretization error');
  plot(h,err(:,1),'r--',h,err(:,2),'b-.',h,err(:,3),'g-',h,err(:,4),'c-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf \eps}');
  ylabel('{\bf Error}');
  legend('L^1-norm A_{up}^{-1}','L^1-norm A_{st}^{-1}','L^1-norm A_{SUPG}^{-1}','Location','NorthWest')
 % p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,1)),1);
%  add_Slope(gca,'SouthEast',p(1),'r--');
 % p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,2)),1);
%  add_Slope(gca,'East',p(1),'b-.');
 % p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,3)),1);
%  add_Slope(gca,'NorthEast',p(1),'g-');
 % p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,4)),1);
%  add_Slope(gca,'Southwest',p(1),'c-');
   
  fig = figure('Name','Discretization error');
  plot(h,err1(:,1),'r--',h,err1(:,2),'b-.',h,err1(:,3),'g-',h,err1(:,4),'c-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf \eps}');
  ylabel('{\bf Error}');
  legend('L^{\infty}-norm A_{up}^{-1}','L^{\infty}-norm A_{st}^{-1}','L^{\infty}-norm A_{SUPG}^{-1}','Location','NorthWest')
  p = polyfit(log(h(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS,1)),1);
  add_Slope(gca,'SouthEast',p(1),'r--');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS,2)),1);
  add_Slope(gca,'East',p(1),'b-.');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS,3)),1);
  add_Slope(gca,'NorthEast',p(1),'g-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS,4)),1);
  add_Slope(gca,'Southwest',p(1),'c-');
  clear all;