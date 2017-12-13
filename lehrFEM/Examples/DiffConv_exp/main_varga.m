% Initialize constants
  JIG=1;
  d1=0.1;                    %impact of supg modification
  d2=0.1;
  a=10^0;                 %amount of diffusivity
  NREFS =5;
  % select test code
  d=getData(1);
     
  QuadRule = P7O6(); 
  
  Mesh.Coordinates =d.Coordinates; 
  Mesh.Elements = d.Elements;
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = d.boundtype;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  err=zeros(NREFS,5);
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
        
   h(i)=get_MeshWidth(Mesh);
   
   %Laplace
   %Al = assemMat_QFE(New_Mesh,@STIMA_Lapl_QFE_lump);
   A = assemMat_QFE(New_Mesh,@STIMA_Lapl_QFE_quad,P3O3());
   %A = assemMat_QFE(New_Mesh,@STIMA_Lapl_QFE);
   weights1=[0 0 0 1/3 1/3 1/3 0];
   [M1,UP1,C1]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights1);
   A_u1=a*A+(M1*UP1);
   L1 = assemLoad_QFE(New_Mesh,P3O3(),d.SOL_Handle,a);
   NCoord=size(Mesh.Coordinates,1);
   NEdg=size(Mesh.Edges,1);
   % Incorporate Dirichlet boundary data
   
   [U_1,FreeDofs] = assemDir_QFE(New_Mesh,[-1 -2],d.U_EX_Handle,a);
   L1 = L1 - A_u1*U_1;
   
   % Solve the linear system
 
   U_1(FreeDofs) = A_u1(FreeDofs,FreeDofs)\L1(FreeDofs);
   
   Vert_Dofs=intersect(FreeDofs,1:NCoord);
   Edge_Dofs=intersect(FreeDofs,NCoord+(1:NEdg));
%    
   M1=A_u1(Vert_Dofs,Vert_Dofs);
   M2=A_u1(Vert_Dofs,Edge_Dofs);
   M3=A_u1(Edge_Dofs,Vert_Dofs);
   M4=A_u1(Edge_Dofs,Edge_Dofs);
   
   W=inv(full(M4))*M3;
   [min(min(-W)) max(max(-W))]
   
   min(min(A_u1*A_u1))
   
   err(i,1) = L2Err_QFE(New_Mesh,U_1,P7O6(),d.U_EX_Handle,0,a);
  
   Dofs(i)=size(New_Mesh.Coordinates,1);
   plot_QFE(U_1,New_Mesh); colorbar
  end; 
  fig = figure('Name','Discretization error');
  plot(h,err(:,1),'ro-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('L^2-error u (up1)','Location','NorthEast')
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,1)),1);
  add_Slope(gca,'SouthEast',p(1),'r-');
  
 clear all;