% Run script discrete differential forms.

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % load test data
    d=getData(1);
    
    % Initialize constant
    EPSI_Handle = @(x,varargin)ones(size(x,1),1);
    MU_HANDLE=@(x,varargin)1;
    ZERO_Handle=@(x,varargin) zeros(size(x,1),2);
    QuadRule=Duffy(TProd(gauleg(0,1,10)));
    NREFS =5;
    JIG =2;  
    
    % Initialize mesh
    Mesh.Coordinates=d.Coordinates;
    Mesh.Elements=d.Elements;
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = d.boundtype;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    % Initialize error-data collecting arrays
    err=zeros(6,NREFS);
    interp=zeros(3,NREFS);
    h=zeros(NREFS,1);
    Dofs=zeros(NREFS,1);
    approx=zeros(NREFS,3);
    
    per = 0;
    progress_bar(per);
    for i = 1:NREFS
        
        if(per < floor(100*i/NREFS))
            per = floor(100*i/NREFS);
            progress_bar(per);  
        end
        
        Mesh=refine_REG(Mesh);
        Mesh=add_Edge2Elem(Mesh);
        Mesh=add_VBdFlags(Mesh);

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
        cVU=assemCochain_0f(New_Mesh,d.V_U_EX_Handle);
        cVUBd=assemBndCochain_0f(New_Mesh,[-2 -3],d.V_U_EX_Handle);
        cCurlU=assemCochain_2f(New_Mesh,d.CURL_U_EX_Handle,QuadRule);
        c1=assemCochain_1f(New_Mesh,d.SOL1_Handle,gauleg(0,1,10));  % -v x curl u righthand side
        c1Bd=assemBndCochain_1f(New_Mesh,[-1], d.SOL1_Handle,gauleg(0,1,10));
        c2=assemCochain_1f(New_Mesh,d.SOL2_Handle,gauleg(0,1,10));  % grad(v.u) righthand side
        c2Bd=assemBndCochain_1f(New_Mesh,[-1],d.SOL2_Handle,gauleg(0,1,10));
        c3=assemCochain_1f(New_Mesh,d.SOL3_Handle,gauleg(0,1,10));  % curl curl u righthand side
        c3Bd=assemBndCochain_1f(New_Mesh, [-1], d.SOL3_Handle,gauleg(0,1,10));
%       L1=assemLoad_W1F(New_Mesh,P7O6(),d.SOL1_Handle);  % -v x curl u righthand side
%       L2=assemLoad_W1F(New_Mesh,P7O6(),d.SOL2_Handle);  % grad(v.u) righthand side
%       L=assemLoad_P0(New_Mesh,P7O6(),d.CURL_U_EX_Handle);  % grad(v.u) righthand side
      
        % contraction
        ContrOne=assemMat_Contr1f(New_Mesh,d.V_Handle);  % contraction of one forms
        V=d.V_Handle(New_Mesh.Coordinates);
        ContrTwo=assemMat_Contr2f(New_Mesh,V);   % contraction of two forms
        
        % toplological derivatives
        TopGrad=assemMat_TopGrad(New_Mesh);   % topological Gradient
        TopRot=assemMat_TopRot(New_Mesh);     % topological Rotation
        
        % several mass matrices
%       M_P0=assemMAT_P0(New_Mesh,@MASS_P0);
%       MassZero=assemMat_Mass0fD(New_Mesh);
       
%       M_LFE=assemMat_LFE(New_Mesh,@MASS_LFE);
        MassOne=assemMat_Mass1fD(New_Mesh);                         % diagonal Mass Matrix of 1-Forms 
        
%       M_W1F = assemMat_W1F(New_Mesh,@MASS_W1F,MU_HANDLE, P3O3());  % Mass Matrix Whitney-1-Forms
%       MassTwo=assemMat_Mass2fD(New_Mesh);                         % diagonal Mass Matrix of Two-Forms      

        % FEM ansatz
        ContrRot=assemMat_W1F(New_Mesh,@STIMA_ContrRot,d.V_Handle, P7O6());  % -v x rot u FEM
        GradContr=assemMat_W1F(New_Mesh,@STIMA_GradContr, d.V_Handle, P7O6());
        
        R=assemMat_W1F(New_Mesh,@STIMA_Curl_W1F,EPSI_Handle,P3O3());
        
        % System Matrix
%       plot_P0(cCurlU,New_Mesh);
%       colorbar;
%       % plot_P0(TopRot*cU,New_Mesh);
%       plot_LFE(ContrOne*cU,New_Mesh); colorbar;
%       plot_LFE((ContrOne*cU+cVUBd),New_Mesh); colorbar;
%       plot_LFE((cVU),New_Mesh); colorbar;
%       plot_Norm_W1F(ContrTwo*TopRot*cU,New_Mesh),colorbar;
%       plot_Norm_W1F(c1,New_Mesh),colorbar;
%       figure;
      % plot_Norm_W1F(TopGrad*ContrOne*cU+TopGrad*cVUBd,New_Mesh);colorbar;
      % plot_Norm_W1F(TopGrad*ContrOne*cU,New_Mesh);colorbar;
      % plot_Norm_W1F(c2,New_Mesh);colorbar;
%       figure;
 
%        c2Bd(FreeDofs)=TopGrad(FreeDofs,:)*ContrOne*cU;
%        plot_Norm_W1F(c2Bd,New_Mesh);colorbar;
%        plot_Norm_W1F(c3,New_Mesh);colorbar; 
        
       % interpolationerror
       %0forms
       interp(1,i)=L2Err_LFE(New_Mesh,cVU,QuadRule,d.V_U_EX_Handle);
       % 1 forms
       interp(2,i)=L2Err_W1F(New_Mesh,cU,QuadRule,d.U_EX_Handle);
       % 2 forms
       interp(3,i)=L2Err_P2f(New_Mesh,cCurlU,QuadRule,d.CURL_U_EX_Handle);
       
       % operator;
       
       err(1,i)=L2Err_W1F(New_Mesh,(TopGrad*ContrOne*cU+TopGrad*cVUBd),QuadRule,d.SOL2_Handle);
       err(2,i)=L2Err_LFE(New_Mesh,ContrOne*cU+cVUBd,QuadRule,d.V_U_EX_Handle);
       err(3,i)=L2Err_W1F(New_Mesh,ContrTwo*TopRot*cU+c1Bd,QuadRule,d.SOL1_Handle);
       err(4,i)=L2Err_W1F(New_Mesh,R*cU,QuadRule,d.SOL3_Handle);
       err(5,i)=L2Err_W1F(New_Mesh,ContrRot*cU,QuadRule, d.SOL1_Handle);
       err(6,i)=L2Err_W1F(New_Mesh,GradContr*cU,QuadRule, d.SOL2_Handle);
       
       Dofs(i)=size(New_Mesh.Edges,1);
        
    end

    fig = figure('Name','Discretization error');
    plot(Dofs,err(1,:),'ro--',Dofs,err(2,:),'go--',Dofs,err(3,:),'bo--',Dofs,err(4,:),'yo--',Dofs,err(5,:),'ko--',Dofs,err(6,:),'co--'); grid('on');
    set(gca,'XScale','log','YScale','log');
    xlabel('{\bf Dofs}');
    ylabel('{\bf Error}');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(1,NREFS-3:NREFS)'),1);
    add_Slope(gca,'NorthEast',p(1),'r-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(2,NREFS-3:NREFS)'),1);
    add_Slope(gca,'East',p(1),'g-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(3,NREFS-3:NREFS)'),1);
    add_Slope(gca,'SouthEast',p(1),'b-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(4,NREFS-3:NREFS)'),1);
    add_Slope(gca,'North',p(1),'y-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(5,NREFS-3:NREFS)'),1);
    add_Slope(gca,'South',p(1),'k-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(6,NREFS-3:NREFS)'),1);
    add_Slope(gca,'NorthWest',p(1),'c-');
    
    %legend('Grad*Conr^1+Grad*Bd','Contr^1+Bd','Conr^2*Rot+Bd','Contr^2+Bd','Grad*Conr^1+Grad*Bd (inn)','','Location','NorthEast') 
    
    legend('Grad*Conr^1+Grad*Bd','Contr^1+Bd','ContrTwo*TopRot','-','Grad*Conr^1+Grad*Bd (inn)','-','Location','NorthEast') 
      
    
    fig = figure('Name','L^2-Interpolation error');
    plot(Dofs,interp(1,:),'ro--',Dofs,interp(2,:),'go--',Dofs,interp(3,:),'bo--'); grid('on');
    set(gca,'XScale','log','YScale','log');
    xlabel('{\bf Dofs}');
    ylabel('{\bf Error}');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(interp(1,NREFS-3:NREFS)'),1);
    add_Slope(gca,'NorthEast',p(1),'r-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(interp(2,NREFS-3:NREFS)'),1);
    add_Slope(gca,'East',p(1),'g-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(interp(3,NREFS-3:NREFS)'),1);
    add_Slope(gca,'SouthEast',p(1),'b-');
    
    legend('0-forms','1-forms','2-forms','Location','NorthEast') 
     
    clear all
    