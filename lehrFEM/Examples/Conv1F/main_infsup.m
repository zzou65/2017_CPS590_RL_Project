% Run script discrete differential forms.

%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    close all;
 %   clear all;
    
    % Initialize constant
    %e_a* (v x curl u)+e_b*(grad v.u)+e_c* curlcurl u +e_m*u
    e_a=1;                            %(v x curl u)   
    e_b=0;                            %(grad v.u)
    e_c=0;                      %curlcurl u 
    e_m=1;                            %u
    %d=getData(12);                   % struct containing Coordinates, Elements and various function handles
    EPSI_Handle = @(x,varargin)ones(size(x,1),1);
    MU_HANDLE = @(x,varargin)1;
    QuadRule = Duffy(TProd(gauleg(0,1,10)));
    l1 = 1; l2 = 0.5;
    V_Handle = @(x,varargin)[l1*ones(size(x,1),1) l2*ones(size(x,1),1)];
    NREFS =4;
    JIG =2;  
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1.5 -1.5; 1.5 -1.5; 1.5 1.5; -1.5 1.5];
    Mesh.Elements = [1 2 3;1 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = [-1 -1 -2 -2];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    err1=zeros(NREFS,1);
    err2=zeros(NREFS,1);
    err3=zeros(NREFS,1);
    err4=zeros(NREFS,1);
    h=zeros(NREFS,1);
    Dofs=zeros(NREFS,1);
    approx=zeros(NREFS,3);
    
    for i = 1:NREFS
    i
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
              
        % contraction
        ContrOne=assemMat_Contr1f(New_Mesh,V_Handle);  % contraction of one forms
        V=V_Handle(New_Mesh.Coordinates);
        ContrTwo=assemMat_Contr2f(New_Mesh,V);   % contraction of two forms
%        ContrRot_UP=assemMat_W1F(New_Mesh,@STIMA_ContrRot_UP,d.V_Handle);  
       
        % topological derivatives
        TopGrad=assemMat_TopGrad(New_Mesh);   % topological Gradient
        TopRot=assemMat_TopRot(New_Mesh);      % topological Rotation
        
    
        % FEM ansatz
        ContrRot=assemMat_W1F(New_Mesh,@STIMA_ContrRot,V_Handle, QuadRule);  % -v x rot u FEM
        GradContr=assemMat_W1F(New_Mesh,@STIMA_GradContr, V_Handle, gauleg(0,1,10));
       
        SUPG=assemMat_W1F(New_Mesh,@STIMA_SUPG_W1F, V_Handle, QuadRule);
        
         % several mass matrices
 
         M_W1F = assemMat_W1F(New_Mesh,@MASS_W1F,MU_HANDLE, P3O3());  % Mass Matrix Whitney-1-Forms
         MassOne=assemMat_Mass1fD(Mesh);                         % diagonal Mass Matrix of 1-Forms 
       
%       M_P0=assemMAT_P0(Mesh,@MASS_P0);
         MassTwo=assemMat_Mass2fD(New_Mesh);                         % diagonal Mass Matrix of Two-Forms    

        
        A=MassOne*ContrTwo*TopRot;              % -v x curl u geom.
        B=MassOne*TopGrad*ContrOne;             % grad(v.u) geom.
        C=TopRot'*MassTwo*TopRot;             % curl curl u geom. 
        
 %       ContrRot_UP=M_W1F*ContrRot_UP;
     
       % standard scheme
       mesh_h=get_MeshWidth(New_Mesh);
       S_st=e_c*C+e_m*M_W1F+e_a*ContrRot+e_b*GradContr;

       % standard scheme
       mesh_h=get_MeshWidth(New_Mesh);
       S_supg=e_c*C+e_m*M_W1F+e_a*ContrRot+e_b*GradContr+mesh_h*SUPG;

       % upwind scheme
       S_up=e_m*M_W1F+e_a*A+e_c*C+e_b*B;
 
       % calculate eigenvalues
       FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
       nonInFlow=find(New_Mesh.BdFlags ~= -1);
       nonOutFlow=find(New_Mesh.BdFlags ~= -2);
       FreeDofs=nonInFlow;
%        E_st = M_W1F(FreeDofs,FreeDofs)\S_st(FreeDofs,FreeDofs); E_st=S_st(FreeDofs,FreeDofs)'*E_st;
%        E_supg = M_W1F(FreeDofs,FreeDofs)\S_supg(FreeDofs,FreeDofs); E_supg=S_supg(FreeDofs,FreeDofs)'*E_supg;
%        E_up = M_W1F(FreeDofs,FreeDofs)\S_up(FreeDofs,FreeDofs); E_up=S_up(FreeDofs,FreeDofs)'*E_up;
       
%        err1(i) = eigs(S_st(FreeDofs,FreeDofs),[],1,'SM');
%        err2(i) = eigs(S_supg(FreeDofs,FreeDofs),[],1,'SM');
%        err3(i) = eigs(S_up(FreeDofs,FreeDofs),[],1,'SM');
       err1(i) = norm(inv(full(S_st(FreeDofs,FreeDofs))));
       err2(i) = norm(inv(full(S_supg(FreeDofs,FreeDofs))));
       err3(i) = norm(inv(full(S_up(FreeDofs,FreeDofs))));

       h(i)=get_MeshWidth(New_Mesh);
       
    end

    fig = figure('Name','Discretization error');
    plot(h,abs(err1),'ro-',h,abs(err2),'go-',h,abs(err3),'bo-'); grid('on');
    set(gca,'XScale','log','YScale','log');
    xlabel('{\bf h}');
    ylabel('{\bf Sangalli-Constant}');
    
    legend('standard','supg','upwind','Location','NorthEast') 
clear all