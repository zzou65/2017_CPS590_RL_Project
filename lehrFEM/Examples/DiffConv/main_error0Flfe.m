% Initialize constants
  
 clear Mesh;
 
  JIG=1;
  d1=1;                    %impact of supg modification
  d2=1;
  a=10^-10                %amount of diffusivity
  NREFS =7;
  % select test code
  d=getData(2);
     
  QuadRule = P7O6(); 
  
  Mesh.Coordinates =d.Coordinates; 
  Mesh.Elements = d.Elements;
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = d.boundtype;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  err=zeros(NREFS,3);
  err1=zeros(NREFS,3);
  h=zeros(NREFS,1);
  Dofs=zeros(NREFS,1);
  
  for i = 1:NREFS
   
   %refine Mesh   
   Mesh = refine_REG(Mesh);
   Mesh=add_Edge2Elem(Mesh);
        
   % Mesh preprocessing
   switch(JIG)
        case 1
          New_Mesh = Mesh;      
        case 2
          Loc = get_BdEdges(Mesh);
          Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
          FixedPos = zeros(size(Mesh.Coordinates,1),1);
          FixedPos(Loc) = 1;
          New_Mesh = jiggle(Mesh,FixedPos);   
        case 3
          Loc = get_BdEdges(Mesh);
          Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
          FixedPos = zeros(size(Mesh.Coordinates,1),1);
          FixedPos(Loc) = 1;
          New_Mesh = smooth(Mesh,FixedPos);
   end
        
   h(i)=get_MeshWidth(New_Mesh);
   
   %Laplace
   A = assemMat_LFE(New_Mesh,@STIMA_Lapl_LFE,P3O3());

   UP4=assemMat_UpLFE(New_Mesh,d.V_Handle);
   B_supg=assemMat_LFE(New_Mesh,@STIMA_SUPG_LFE,P3O3(), d.V_Handle,a,d1,d2);
   B=assemMat_LFE(New_Mesh,@STIMA_Conv_LFE,d.V_Handle,P7O4());

   A_u4=a*A+(UP4);
   A_supg=a*A+B+B_supg;
   A_s=a*A+B;
   L = assemLoad_LFE(New_Mesh,P3O3(),d.SOL_Handle,a);
   L_supg=assemLoad_LFE_SUPG(New_Mesh,P3O3(),d.V_Handle,d.SOL_Handle,a,d1,d2);
   
   %Incorporate Dirichlet boundary data
   
   [U_s,FreeDofs] = assemDir_LFE(New_Mesh,[-1 -2],d.U_EX_Handle,a);
   U_u4=U_s;
   U_supg=U_s;
    
   L_u4 = L - A_u4*U_u4;
   L_s = L - A_s*U_s;
   L_supg = L+L_supg - A_supg*U_supg;
   
   % Solve the linear system
 
   U_u4(FreeDofs) = A_u4(FreeDofs,FreeDofs)\L_u4(FreeDofs);
   U_s(FreeDofs) = A_s(FreeDofs,FreeDofs)\L_s(FreeDofs);
   U_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);
   %
   
   err1(i,1) = L2Err_LFE(New_Mesh,U_u4,P7O6(),d.U_EX_Handle,0,a);
   err1(i,2) = L2Err_LFE(New_Mesh,U_s,P7O6(),d.U_EX_Handle,0,a);
   err1(i,3) = L2Err_LFE(New_Mesh,U_supg,P7O6(),d.U_EX_Handle,0,a)
   %
   Dofs(i)=size(New_Mesh.Coordinates,1); 
  % plot_LFE(U_u4,New_Mesh); colorbar
   % 
 end; 
 
 fig = figure('Name','Discretization error');
 plot(h,err1(:,1),'go-',h,err1(:,2),'bo-',h,err1(:,3),'yo-'); grid('on');
 set(gca,'XScale','log','YScale','log');
 xlabel('{\bf h}');
 ylabel('{\bf Error}');
 legend('L^2-error u (up1)','L^2-error u','L^2-error u(supg)','Location','NorthEast');
 p = polyfit(log(h(1:NREFS)),log(err1(1:NREFS,1)),1);
 add_Slope(gca,'West',p(1),'g-');
 p = polyfit(log(h(1:NREFS)),log(err1(1:NREFS,2)),1);
 add_Slope(gca,'East',p(1),'b-');
 p = polyfit(log(h(1:NREFS)),log(err1(1:NREFS,3)),1);
 add_Slope(gca,'Southwest',p(1),'y-');
 
 % print -depsc converge2_comp2_-10.eps
 %   
%   fig = figure('Name','Discretization error');
%   plot(h,err(:,4),'go-',h,err(:,1),'ro-',h,err(:,2),'bo-',h,err(:,3),'yo-',h,err(:,5),'co-',h,err(:,6),'mo-'); grid('on');
%   set(gca,'XScale','log','YScale','log');
%   xlabel('{\bf h}');
%   ylabel('{\bf Error}');
%   legend('L^{\infty}-error u (up1)','L^{\infty}-error u (up2)','L^{\infty}-error u (up3)',...
%       'L^{\infty}-error u (up4)','L^{\infty}-error u','L^{\infty}-error u(supg)','Location','NorthEast');
%   p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,4)),1);
%   add_Slope(gca,'West',p(1),'g-');
%   p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,1)),1);
%   add_Slope(gca,'SouthEast',p(1),'r-');
%   p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,2)),1);
%   add_Slope(gca,'East',p(1),'b-');
%   p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,3)),1);
%   add_Slope(gca,'Southwest',p(1),'y-');
%   p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,5)),1);
%   add_Slope(gca,'NorthWest',p(1),'c-');
%   p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,6)),1);
%   add_Slope(gca,'NorthEast',p(1),'m-');
 % print -depsc converge2_part2_-10.eps
  
%clear all;