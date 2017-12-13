% Run script discrete differential forms.

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % Initialize constant
    EPSI_Handle = @(x,varargin)ones(size(x,1),1);
    MU_HANDLE=@(x,varargin)1;
    Zero_Handle=@(x,varargin) zeros(size(x,1),2);
        
%     % Data I
%     
%      V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)];
%     
%      U_EX_Handle=@(x,varargin)[(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1)) ...
%                                (1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2))];
%      SOL1_Handle=@(x,varargin)[  -x(:,2).*sin(pi.*x(:,1))+x(:,1).*sin(pi.*x(:,2)) ...
%                               +2.*x(:,2).*sin(pi.*x(:,1))-2.*x(:,1).*sin(pi.*x(:,2))];    % - v x curl u
%      SOL2_Handle=@(x,varargin)[-pi.*(-1+x(:,2).^2).*cos(pi.*x(:,1))-x(:,1).*sin(pi.*x(:,2)) ...
%                                -pi/2.*(-1+x(:,1).^2).*cos(pi.*x(:,2))-2*x(:,2).*sin(pi.*x(:,1))]; % grad (v u)
%      SOL3_Handle=@(x,varargin)[2.*(-pi.*x(:,1).*cos(pi.*x(:,2))+sin(pi.*x(:,1))) ...
%                                2.*(-pi.*x(:,2).*cos(pi.*x(:,1))+sin(pi.*x(:,2)))]; % curl curl u
%      CURL_U_EX_Handle=@(x,varargin)2.*x(:,2).*sin(pi.*x(:,1))- 2.*x(:,1).*sin(pi.*x(:,2));    % curl u

%   Data II


    V_Handle=@(x,varargin)ones(size(x,1),2);
% 
%     U_EX_Handle=@(x,varargin)[sin(pi.*x(:,2)).*sin(pi.*x(:,1)) ...
%                               sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
%     CURL_U_EX_Handle=@(x,varargin)[-pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))+pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2))]  ;                    
%     SOL1_Handle=@(x,varargin)[+pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
%                               -pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))+pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2))];    % - v x curl u
%     SOL2_Handle=@(x,varargin)[2.*pi.*sin(pi.*x(:,2)).*cos(pi.*x(:,1))...
%                               2.*pi.*sin(pi.*x(:,1)).*cos(pi.*x(:,2))]; % grad (v u)
%     SOL3_Handle=@(x,varargin)[pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
%                              +pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2))]; % curl curl u
%     
  %  Data III

%     U_EX_Handle=@(x,varargin)[1-x(:,1).^2 1-x(:,2).^2];
%     SOL1_Handle=@(x,varargin) zeros(size(x,1),2);
%     SOL2_Handle=@(x,varargin) [-2.*x(:,1) -2.*x(:,2)];
%     SOL3_Handle=@(x,varargin) zeros(size(x,1), 2);
%     CURL_U_EX_Handle=@(x,varargin)zeros(size(x,1),1)  ;
    
 %  Data IV
% 
    U_EX_Handle=@(x,varargin)[1-x(:,1) 1-x(:,2)];
    SOL1_Handle=@(x,varargin) zeros(size(x,1),2);
    SOL2_Handle=@(x,varargin) -ones(size(x(:,1),1),2);
    SOL3_Handle=@(x,varargin) zeros(size(x,1), 2);
    CURL_U_EX_Handle=@(x,varargin)zeros(size(x,1),1)  ;    
                              
           

    QuadRule=Duffy(TProd(gauleg(0,1,4)));
    NREFS =6;
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 4;2 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = [-1 -1 -2 -2];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    err=zeros(NREFS,1);
    err2=zeros(NREFS,1);
    h=zeros(NREFS,1);
    Dofs=zeros(NREFS,1);
    approx=zeros(NREFS,3);
    
    for i = 1:NREFS
       
        Mesh=refine_REG(Mesh);
        Mesh = add_Edge2Elem(Mesh);
        
        L1=assemLoad_W1F(Mesh,P7O6(),SOL1_Handle);  % -v x curl u righthand side
        L2=assemLoad_W1F(Mesh,P7O6(),SOL2_Handle);  % grad(v.u) righthand side
        L3=assemLoad_W1F(Mesh,P7O6(),SOL3_Handle);  % curl curl u righthand side
        L4=assemLoad_W1F(Mesh,P7O6(),U_EX_Handle);  % u righthand side
        
        ContrOne=assemMat_ContrOne(Mesh,V_Handle);   % contraction of one forms
        ContrTwo=assemMat_ContrTwo(Mesh,V_Handle);   % contraction of two forms
       
        TopGrad=assemMat_TopGrad(Mesh);   % topological Gradient
        TopRot=assemMat_TopRot(Mesh);     % topological Rotation
        
        MassOne=assemMat_MassOneD(Mesh);                     % diagonal ...Mass Matrix of zero-Forms 
        M = assemMat_W1F(Mesh,@MASS_W1F,MU_HANDLE, P7O6());  % Mass Matrix Whitney-1-Forms
        MassTwo=assemMat_MassTwoD(Mesh);                     % diagonal Mass Matrix of Two-Forms      
       
        A=M*ContrTwo*TopRot;              % -v x curl u geom.
        B=M*TopGrad*ContrOne;             % grad(v.u) geom.
        C=TopRot'*MassTwo*TopRot;         % curl curl u geom. 
%       R=assemMat_W1F(Mesh,@STIMA_Curl_W1F,EPSI_Handle,P3O3())  % curl curl u FEM;

        % System Matrix
        
        S =(A+C+M);
        
        %  Direchlet Boundary data
        
        [U,FreeDofs] = assemDir_W1F(Mesh,[-1 -2],U_EX_Handle,gauleg(0,1,4));
   
        L=L1+L3+L4;
        
        L = L - S*U;
        
        % Solve the system
     
        U(FreeDofs) = S(FreeDofs,FreeDofs)\L(FreeDofs);
        plot_Norm_W1f(U,Mesh);
        err(i)=L2Err_W1F(Mesh,U,QuadRule,U_EX_Handle)
        err2(i)=HCurlSErr_W1F(Mesh,-U,QuadRule,CURL_U_EX_Handle)
        h(i)=get_MeshWidth(Mesh);
        Dofs(i)=size(Mesh.Edges,1);
       
        end

    fig = figure('Name','Discretization error');
    plot(Dofs,err,'ro--',Dofs,err2,'go--'); grid('on');
    set(gca,'XScale','log','YScale','log');
    xlabel('{\bf Dofs}');
    ylabel('{\bf Error}');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(NREFS-3:NREFS)),1);
    add_Slope(gca,'East',p(1),'r-');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS)),1);
    add_Slope(gca,'SouthEast',p(1),'g-');
    
    legend('L^2','Hcurl-semi','Location','NorthEast') 
 clear all
    