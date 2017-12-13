% Run script discrete differential forms.

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % Initialize constant
    EPSI_Handle = @(x,varargin)ones(size(x,1),1);
    MU_HANDLE=@(x,varargin)1;
    ZERO_Handle=@(x,varargin) zeros(size(x,1),2);
        
    QuadRule=Duffy(TProd(gauleg(0,1,10)));
    NREFS =5;
    JIG = 3;  
    d=getData(1);
        
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
    h=zeros(NREFS,1);
    Dofs=zeros(NREFS,1);
    approx=zeros(NREFS,3);
    
    for i = 1:NREFS
    
        Mesh=refine_REG(Mesh);
        Mesh = add_Edge2Elem(Mesh);
        
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
%       cVU=assemCochain_0f(New_Mesh,d.V_U_EX_Handle);
%       cVUBd=assemBndCochain_0f(New_Mesh,[-2 -3],d.V_U_EX_Handle);
%       cCurlU=assemCochain_2f(New_Mesh,d.CURL_U_EX_Handle,QuadRule);
%       c1=assemCochain_1f(New_Mesh,d.SOL1_Handle,gauleg(0,1,10));  % -v x curl u righthand side
%       c1Bd=assemBndCochain_1f(New_Mesh,[-1], d.SOL1_Handle,gauleg(0,1,10));
%       c2=assemCochain_1f(New_Mesh,d.SOL2_Handle,gauleg(0,1,10));  % grad(v.u) righthand side
%       c2Bd=assemBndCochain_1f(New_Mesh,[-1],d.SOL2_Handle,gauleg(0,1,10));
%       c3=assemCochain_1f(New_Mesh,d.SOL3_Handle,gauleg(0,1,10));  % curl curl u righthand side
%       c3Bd=assemBndCochain_1f(New_Mesh, [-1], d.SOL3_Handle,gauleg(0,1,10));
      
        % contraction
        ContrOne=assemMat_Contr1f(New_Mesh,d.V_Handle);  % contraction of one forms
        V=d.V_Handle(New_Mesh.Coordinates);
        ContrTwo=assemMat_Contr2f(New_Mesh,V);   % contraction of two forms
       
        % topological derivatives
        TopGrad=assemMat_TopGrad(New_Mesh);   % topological Gradient
        TopRot=assemMat_TopRot(New_Mesh);     % topological Rotation
        
         % several mass matrices
%       M_LFE=assemMat_LFE(Mesh,@MASS_LFE);
         MassZero=assemMat_Mass0fD(Mesh);
 
         M_W1F = assemMat_W1F(New_Mesh,@MASS_W1F,MU_HANDLE, P3O3());  % Mass Matrix Whitney-1-Forms
         MassOne=assemMat_Mass1fD(Mesh);                         % diagonal Mass Matrix of 1-Forms 
       
         MassTwo=assemMat_Mass2fD(New_Mesh);                         % diagonal Mass Matrix of Two-Forms    

        
        A = M_W1F*ContrTwo*TopRot;              % -v x curl u geom.
        B = M_W1F*TopGrad*ContrOne;             % grad(v.u) geom.
        C =TopRot'*MassTwo*TopRot;              % curl curl u geom. 
 
        % System Matrix
        
        S=A+B;
                
        % boundary data on inflow boundary and righthandside
        
        FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
        nonInFlow=find(New_Mesh.BdFlags ~=-1);
        nonOutFlow=find(New_Mesh.BdFlags ~=-2);
        
        c1=assemCochain_1f(New_Mesh,d.SOL1_Handle,gauleg(0,1,10));  % -v x curl u righthand side
        c2=assemCochain_1f(New_Mesh,d.SOL2_Handle,gauleg(0,1,10));  % grad(v.u) righthand side
        c3=assemCochain_1f(New_Mesh,d.SOL3_Handle,gauleg(0,1,10));  % curl curl rightand side
        c4=assemCochain_1f(New_Mesh,d.U_EX_Handle,gauleg(0,1,10));  % u righthsnd side        
        cVUBd=assemBndCochain_0f(New_Mesh,[-2 -3],d.V_U_EX_Handle);
         
        cU=assemCochain_1f(New_Mesh,d.U_EX_Handle,gauleg(0,1,10));
        U=cU;
        %FreeDofs=nonInFlow;
        U(FreeDofs)=0;
        cf=c1+c2;
        h1=c1;
        h1(nonInFlow)=0;
      
        %Lf=M_W1F*(cf-h1)-S*U;
        %Lf=M_W1F*(cf-TopGrad*cVUBd)-S*U;
        Lf=M_W1F*(cf-TopGrad*cVUBd-h1)-S*U;
        
        %Lf=M_W1F*(cf-h)-(C+M_W1F)*M_W1F*U;
        %Lf=M_W1F*(cf)-S*U;
        U(FreeDofs)=S(FreeDofs,FreeDofs)\Lf(FreeDofs);
        
        plot_Norm_W1F(U,New_Mesh); colorbar; 
        plot_Norm_W1F(cU,New_Mesh); colorbar;
        err1(i)=L2Err_W1F(New_Mesh,U,QuadRule,d.U_EX_Handle)
        %err2(i)=L2Err_W1Finn(New_Mesh,U,QuadRule,U_EX_Handle)
        err2(i)=HCurlSErr_W1F(New_Mesh,U,QuadRule,d.CURL_U_EX_Handle)
        h(i)=get_MeshWidth(New_Mesh);
        Dofs(i)=size(New_Mesh.Edges,1);
   end

    fig = figure('Name','Discretization error');
    plot(Dofs,err1,'ro--',Dofs,err2,'go--'); grid('on');
    set(gca,'XScale','log','YScale','log');
    xlabel('{\bf Dofs}');
    ylabel('{\bf Error}');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS)),1);
    add_Slope(gca,'East',p(1),'r-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS)),1);
    add_Slope(gca,'SouthEast',p(1),'g-');
    
    legend('L^2','Hcurl-semi','Location','NorthEast') 
% clear all