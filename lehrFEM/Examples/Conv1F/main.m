% Run script discrete differential forms.

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    %close all;
    %   clear all;
    
    % Initialize constant
    %e_a* (v x curl u)+e_b*(grad v.u)+e_c* curlcurl u +e_m*u
    e_a=0;                            % (v x curl u)   
    e_b=0;                            % (grad v.u)
    e_c=1;                     % curlcurl u 
    e_m=1;                           %u
    e_s=0;                            %supg
    d=getData(1);                % struct containing Coordinates, Elements and various function handles
    EPSI_Handle = @(x,varargin)ones(size(x,1),1);
    MU_HANDLE=@(x,varargin)1;
    QuadRule=Duffy(TProd(gauleg(0,1,20)));
    NREFS =5;
    JIG = 2;  
    
    % Initialize mesh
    
    Mesh.Coordinates = d.Coordinates;
    Mesh.Elements = d.Elements;
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = d.boundtype;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    err1=zeros(NREFS,1);
    err2=zeros(NREFS,1);
    err3=zeros(NREFS,1);
    err4=zeros(NREFS,1);
    h=zeros(NREFS,1);
    Dofs=zeros(NREFS,1);
    approx=zeros(NREFS,3);
    
    for i = 1:NREFS
    
        Mesh=refine_REG(Mesh);
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
        
        cU=assemCochain_1f(New_Mesh,d.U_EX_Handle,gauleg(0,1,10));
%      cVU=assemCochain_0f(New_Mesh,d.V_U_EX_Handle);
%      cVUBd=assemBndCochain_0f(New_Mesh,[-2 -3],d.V_U_EX_Handle);
%       cCurlU=assemCochain_2f(New_Mesh,d.CURL_U_EX_Handle,QuadRule);
%       c1=assemCochain_1f(New_Mesh,d.SOL1_Handle,gauleg(0,1,10));  % -v x curl u righthand side
%       c1Bd=assemBndCochain_1f(New_Mesh,[-1], d.SOL1_Handle,gauleg(0,1,10));
        c2=assemCochain_1f(New_Mesh,d.SOL2_Handle,gauleg(0,1,10));  % grad(v.u) righthand side
%       c2Bd=assemBndCochain_1f(New_Mesh,[-1],d.SOL2_Handle,gauleg(0,1,10));
%       c3=assemCochain_1f(New_Mesh,d.SOL3_Handle,gauleg(0,1,10));  % curl curl u righthand side
%       c3Bd=assemBndCochain_1f(New_Mesh, [-1], d.SOL3_Handle,gauleg(0,1,10));
      
        % contraction
        ContrOne=assemMat_Contr1f(New_Mesh,d.V_Handle);  % contraction of one forms
        V=d.V_Handle(New_Mesh.Coordinates);
        ContrTwo=assemMat_Contr2f(New_Mesh,V);   % contraction of two forms
        ContrRot_UP=assemMat_W1F(New_Mesh,@STIMA_ContrRot_UP,d.V_Handle);  
       
        % topological derivatives
        TopGrad=assemMat_TopGrad(New_Mesh);   % topological Gradient
        TopRot=assemMat_TopRot(New_Mesh);      % topological Rotation
        
        % direct LieDerivative
        Lie=assemMat_W1F(New_Mesh,@STIMA_Lie_W1F, d.V_Handle);
        
        % FEM ansatz
        ContrRot=assemMat_W1F(New_Mesh,@STIMA_ContrRot,d.V_Handle, QuadRule);  % -v x rot u FEM
        ContrRot2=assemMat_W1F(Mesh,@STIMA_ContrRot,d.V_Handle, QuadRule);  % -v x rot u FEM
        GradContr=assemMat_W1F(New_Mesh,@STIMA_GradContr, d.V_Handle, gauleg(0,1,10));
       % GradContr=assemMat_W1F(Mesh,@STIMA_GradContr, d.V_Handle, gauleg(0,1,10));
       
        
        %GradContr2=assemMat_W1F(New_Mesh,@STIMA_GradContr2, d.V_Handle, QuadRule);
        SUPG=assemMat_W1F(New_Mesh,@STIMA_SUPG_W1F, d.V_Handle, QuadRule);
        
         % several mass matrices
         M_LFE=assemMat_LFE(Mesh,@MASS_LFE);
         MassZero=assemMat_Mass0fD(Mesh); 
  
         M_W1F = assemMat_W1F(New_Mesh,@MASS_W1F,MU_HANDLE, P3O3());  % Mass Matrix Whitney-1-Forms
         MassOne=assemMat_Mass1fD(Mesh);                         % diagonal Mass Matrix of 1-Forms 
       
%       M_P0=assemMAT_P0(Mesh,@MASS_P0);
         MassTwo=assemMat_Mass2fD(New_Mesh);                         % diagonal Mass Matrix of Two-Forms    

        
        A=M_W1F*ContrTwo*TopRot;              % -v x curl u geom.
        B=M_W1F*TopGrad*ContrOne;             % grad(v.u) geom.
        C=TopRot'*MassTwo*TopRot;             % curl curl u geom. 
        
        ContrRot_UP=ContrRot_UP;
 
        % boundary data on inflow boundary and righthandside
        
        c1=assemCochain_1f(New_Mesh,d.SOL1_Handle,gauleg(0,1,10));     % -v x curl u righthand side
        c2=assemCochain_1f(New_Mesh,d.SOL2_Handle,gauleg(0,1,10));     % grad(v.u) righthand side
        c3=assemCochain_1f(New_Mesh,d.SOL3_Handle,gauleg(0,1,10));     % curl curl rightand side
        c4=assemCochain_1f(New_Mesh,d.U_EX_Handle,gauleg(0,1,10));      % u righthand side        
        cVUBd=assemBndCochain_0f(New_Mesh,[-2 -3],d.V_U_EX_Handle); % boundary data on inflow boundary
        c1Bd=assemBndCochain_1f(New_Mesh,[-2],d.SOL1_Handle,gauleg(0,1,10)); % boundary 

        % Stabelization
       % InvM0=inv(MassZero);
       % Stab=M_W1F*TopGrad*InvM0*TopGrad'*M_W1F;
        
        % exact solution
        cU=assemCochain_1f(New_Mesh,d.U_EX_Handle,gauleg(0,1,10));
%         plot_Norm_W1F(cU,New_Mesh); colorbar;
%         
        % calculate solution

%         % System Matrix
%         S=M_W1F;
%         FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
%         nonInFlow=find(New_Mesh.BdFlags ~=-1);
%         cf=c4;
%         cU(FreeDofs)=0;
%         Lf=M_W1F*cf-S*cU;
        
%         % System Matrix
%         S=M_W1F+C;
%         FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
%         nonInFlow=find(New_Mesh.BdFlags ~=-1);
%         cf=c4+c3;
%         cU(FreeDofs)=0;
%         Lf=M_W1F*cf-S*cU;
        
%         % System Matrix
%         S=B+C;
%         FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
%         nonInFlow=find(New_Mesh.BdFlags ~=-1);
%         %FreeDofs=nonInFlow;
%         cf=c2+c3;
%         cU(FreeDofs)=0;
%         Lf=M_W1F*(cf-TopGrad*cVUBd)-S*cU;
        
        
%         % System Matrix
%         S=B+A;
%         nonInFlow=find(New_Mesh.BdFlags ~=-1);
%         FreeDofs=nonInFlow;
%         cf=c2+c1;
%         cU(FreeDofs)=0;
%         Lf=M_W1F*(cf-c1Bd-TopGrad*cVUBd)-S*cU;
%         

%        % System Matrix
%        S=M_W1F+ContrRot;
%        FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
%        cf=c4+c1;
%        cU(FreeDofs)=0;
%        Lf=M_W1F*cf-S*cU;

%        % System Matrix
%        S=M_W1F+GradContr;
%        FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
%        cf=c4+c2;
%        cU(FreeDofs)=0;
%        Lf=M_W1F*cf-S*cU;

       % System Matrix
%        cU_st=cU;
%        S_st=e_c*C+e_m*M_W1F+ContrRot;
%        FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
% %        nonInFlow=find(New_Mesh.BdFlags ~=-1);
% %        FreeDofs=nonInFlow;
%        cf=e_c*c3+e_m*c4+c1;
%        cU_st(FreeDofs)=0;
%        Lf_st=M_W1F*cf-S_st*cU_st;
        
%        % System Matrix
%        cU_up=cU;
%        S_up=e_c*C+e_m*M_W1F+A;
%        FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
%        cf=e_c*c3+e_m*c4+c1;
%        cU_up(FreeDofs)=0;
%        Lf_up=M_W1F*(cf-c1Bd)-S_up*cU_up;
           
       % System Matrix and Righthand side
       nonInFlow=find(New_Mesh.BdFlags ~= -1);
       nonOutFlow=find(New_Mesh.BdFlags ~= -2);
       
       % standard scheme
       mesh_h=get_MeshWidth(New_Mesh);
       cU_st=cU;
       S_st=e_c*(C)+e_m*M_W1F+e_a*ContrRot+e_b*GradContr+e_s*mesh_h*SUPG;
       cf=e_c*c3+e_a*c1+e_m*c4+e_b*c2;
       
       FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
     %  FreeDofs=nonInFlow;
       
       cU_st(FreeDofs)=0;
       Lf_st=M_W1F*(cf)-S_st*cU_st;
       
       cU_st(FreeDofs)=S_st(FreeDofs,FreeDofs)\Lf_st(FreeDofs);
       
       % upwind scheme
       nonInFlow=find(New_Mesh.BdFlags ~= -1);
       nonOutFlow=find(New_Mesh.BdFlags ~= -2);
    
       cU_up=cU;
      % S_up=e_m*speye(size(Mesh.Edges(),1))+e_a*ContrRot_UP+e_c*C+e_b*B;
       S_up=e_m*M_W1F+e_a*A+e_c*(C)+e_b*B;
       cf=e_m*c4+e_a*c1+e_c*c3+e_b*c2;
       h1=c1;
        
     %FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
       FreeDofs=nonInFlow;
       h1(nonInFlow)=0;
       %FreeDofs=nonOutFlow;
       cU_up(FreeDofs)=0;
       Lf_up=M_W1F*(cf-e_a*h1)-S_up*cU_up;
       %Lf_up=M_W1F*(cf-e_a*c1Bd-e_b*TopGrad*cVUBd)-S_up*cU_up;

       % solve equations
     
       cU_up(FreeDofs)=cU_up(FreeDofs)+S_up(FreeDofs,FreeDofs)\Lf_up(FreeDofs);
 
       plot_Norm_W1F(cU,New_Mesh);
       plot_Norm_W1F(cU_st,New_Mesh);
       plot_Norm_W1F(cU_up,New_Mesh);
      %  plot_LFE(ContrOne*(Lf_up),New_Mesh);  colorbar;
       plot_LFE(ContrOne*(cU_up),New_Mesh); colorbar;
     %  plot_LFE(ContrOne*(cU_st),New_Mesh); colorbar;
       
       err1(i)=L2Err_W1F_mod(New_Mesh,cU_st,QuadRule,d.U_EX_Handle);
       err2(i)=HCurlSErr_W1F_mod(New_Mesh,cU_st,QuadRule,d.CURL_U_EX_Handle);
       err3(i)=L2Err_W1F_mod(New_Mesh,cU_up,QuadRule,d.U_EX_Handle);
       err4(i)=HCurlSErr_W1F_mod(New_Mesh,cU_up,QuadRule,d.CURL_U_EX_Handle);

       h(i)=get_MeshWidth(New_Mesh);
       Dofs(i)=size(New_Mesh.Edges,1);
       
    end

    fig = figure('Name','Discretization error');
    plot(Dofs,err1,'ro--',Dofs,err2,'go--',Dofs,err3,'bo--',Dofs,err4,'yo--'); grid('on');
    set(gca,'XScale','log','YScale','log');
    xlabel('{\bf Dofs}');
    ylabel('{\bf Error}');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS)),1);
    add_Slope(gca,'SouthWest',p(1),'r-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS)),1);
    add_Slope(gca,'SouthEast',p(1),'g-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err3(NREFS-3:NREFS)),1);
    add_Slope(gca,'SouthWest',p(1),'b-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err4(NREFS-3:NREFS)),1);
    add_Slope(gca,'SouthEast',p(1),'y-');
    
    legend('L^2','Hcurl-semi','L^2(up)','Hcurl-semi(up)','Location','NorthEast') 
clear all