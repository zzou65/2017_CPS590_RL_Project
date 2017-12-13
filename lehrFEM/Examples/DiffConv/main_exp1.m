% solves a convection diffusion equation for different sets of parameters

%   Copyright 2008-2008 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  clear all
 clear Mesh;
 
  JIG=2;
  d1=1;                    %impact of supg modification
  d2=1;
  NREFS =4;                % number of refinement steps
  aa = 10.^-(0:2:10);      % different diffusivity constants
  NPREREFS = 1;            % number of pre-refinements
  
    
  
  % select test code
  d=getData(3);
     
  QuadRule = P7O6(); 
  
  Mesh.Coordinates =d.Coordinates; 
  Mesh.Elements = d.Elements;
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = d.boundtype;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  for i=1:NPREREFS
      Mesh = refine_REG(Mesh);
  end
  
  exp1_errStGal=zeros(NREFS,length(aa));
  exp1_errUPQuad=zeros(NREFS,length(aa));
  exp1_errSUPG=zeros(NREFS,length(aa));
  exp1_errFVol=zeros(NREFS,length(aa));

  exp1_h=zeros(NREFS,1);
  exp1_Dofs=zeros(NREFS,1);
  
  exp1_diffrange = aa;
  
  
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
   New_Mesh = add_MidPoints(New_Mesh, 'barycentric');
        
   exp1_h(i)=get_MeshWidth(New_Mesh);
   exp1_Dofs(i)=size(New_Mesh.Coordinates,1); 
   
   exp1_meshes(i) = New_Mesh;
   
   %Laplace
   
   A = assemMat_LFE(New_Mesh,@STIMA_Lapl_LFE,P3O3());

   UP4=assemMat_upLFEVert(New_Mesh,d.V_Handle);

   B=assemMat_LFE(New_Mesh,@STIMA_Conv_LFE,d.V_Handle,P7O4());
   
   j=0;
   for a = aa
       j=j+1;
       disp(['Running test ' num2str(j+length(aa)*(i-1)) ' of ' num2str(length(aa)*NREFS)]);
       
       B_supg=assemMat_LFE(New_Mesh,@STIMA_SUPG_LFE,P3O3(), d.V_Handle,a,d1,d2);

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
       
       % finite volume discretization
       
       k = @(x,varargin)a;
       c = @(x,varargin)[1 0];
       r = [];
       df = @(x,varargin)x(:,1)-1./(1-exp(-1/a)).*(exp((x(:,1)-1)/a)-exp(-1/a));
       n = [];
       f = @(x,varargin)1;


       [A_fvol, U_fvol, L_fvol, fd_fvol] = assemMat_CD_LFV(New_Mesh, k, c, r, df, n, f);



       % Solve the linear systems

       U_u4(FreeDofs) = A_u4(FreeDofs,FreeDofs)\L_u4(FreeDofs);
       U_s(FreeDofs) = A_s(FreeDofs,FreeDofs)\L_s(FreeDofs);
       U_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);
       U_fvol(fd_fvol) = A_fvol(fd_fvol,fd_fvol)\L_fvol(fd_fvol);
       
       exp1_solStGal(i,j).sol=U_s;
       exp1_solUP(i,j).sol=U_u4;
       exp1_solSUPG(i,j).sol=U_supg;
       exp1_solFVol(i,j).sol=U_fvol;
       
       % compute errors

       exp1_errStGal(i,j) =  L2Err_LFE(New_Mesh,U_s,P7O6(),d.U_EX_Handle,0,a);
       exp1_errUPQuad(i,j) = L2Err_LFE(New_Mesh,U_u4,P7O6(),d.U_EX_Handle,0,a);
       exp1_errSUPG(i,j) =   L2Err_LFE(New_Mesh,U_supg,P7O6(),d.U_EX_Handle,0,a);
       exp1_errFVol(i,j) =   L2Err_LFE(New_Mesh,U_fvol, P7O6(), d.U_EX_Handle,0,a);
   end 
 end
 
 save 'experiment1.mat' exp1_*

  
clear all;
