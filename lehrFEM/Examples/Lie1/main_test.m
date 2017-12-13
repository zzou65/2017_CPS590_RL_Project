% Run script discrete differential forms.

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % Initialize constant
    EPSI_Handle = @(x,varargin)ones(size(x,1),1);
    SIGMA_Handle = @(x,varargin)ones(size(x,1),1);
    MU_HANDLE=@(x,varargin)1;
    Zero_Handle=@(x,varargin) zeros(size(x,1),2);
        
    V_Handle=@(x,varargin)ones(size(x,1),2);
    Grad_V_Handle=@(x,varargin)zeros(size(x,1),4);

    % Data I
    
%      V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)];
%     
%      U_EX_Handle=@(x,varargin)[(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1)) ...
%                                (1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2))];
%      V_U_EX_Handle=@(x,varargin)(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1))+0.5.*(1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2));                      
%      SOL1_Handle=@(x,varargin)[  -x(:,2).*sin(pi.*x(:,1))+   x(:,1).*sin(pi.*x(:,2)) ...
%                               +2.*x(:,2).*sin(pi.*x(:,1))-2.*x(:,1).*sin(pi.*x(:,2))];    % - v x curl u
%      SOL2_Handle=@(x,varargin)[-pi.*(-1+x(:,2).^2).*cos(pi.*x(:,1))-x(:,1).*sin(pi.*x(:,2)) ...
%                                -pi/2.*(-1+x(:,1).^2).*cos(pi.*x(:,2))-2*x(:,2).*sin(pi.*x(:,1))]; % grad (v u)
%      SOL3_Handle=@(x,varargin)[2.*(-pi.*x(:,1).*cos(pi.*x(:,2))+sin(pi.*x(:,1))) ...
%                                2.*(-pi.*x(:,2).*cos(pi.*x(:,1))+sin(pi.*x(:,2)))]; % curl curl u

    % Data II
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
 %  Data III

%    U_EX_Handle=@(x,varargin)[1-x(:,1).^2 1-x(:,2).^2];
%    SOL1_Handle=@(x,varargin) zeros(size(x,1),2);
%    SOL2_Handle=@(x,varargin) [-2.*x(:,1) -2.*x(:,2)];
%    SOL3_Handle=@(x,varargin) zeros(size(x,1), 2);
%    CURL_U_EX_Handle=@(x,varargin)zeros(size(x,1),1)  ;
%    UU_EX_Handle=@(x,varargin)[(1-x(:,1).^2).^2+( 1-x(:,2).^2).^2];
%     U_EX_Handle=@(x,varargin) 1/sqrt(2).*ones(size(x,1),2);
%     UU_EX_Handle=@(x,varargin)ones(size(x,1),1);
%     
 %  Data IV

%    V_Handle=@(x,varargin)[ones(size(x,1),1) 1*ones(size(x,1),1)];
%    U_EX_Handle=@(x,varargin)[1-x(:,1) 1-x(:,2)];
%    SOL1_Handle=@(x,varargin) zeros(size(x,1),2);
%    SOL2_Handle=@(x,varargin) [-ones(size(x(:,1),1),1) -1*ones(size(x(:,1),1),1)];
%    SOL3_Handle=@(x,varargin) zeros(size(x,1), 2);
%    CURL_U_EX_Handle=@(x,varargin)zeros(size(x,1),1)  ;    

%  Data V

    V_Handle=@(x,varargin)[ones(size(x,1),1) 1*ones(size(x,1),1)];
    U_EX_Handle=@(x,varargin)[-x(:,2) +x(:,1)];
    CURL_U_EX_Handle=@(x,varargin)2.*ones(size(x,1),1);
    SOL1_Handle=@(x,varargin) [-2*ones(size(x(:,1),1),1) 2*ones(size(x(:,1),1),1)]
    SOL2_Handle=@(x,varargin) [2*ones(size(x(:,1),1),1) 2*ones(size(x(:,1),1),1)];
    SOL3_Handle=@(x,varargin) zeros(size(x,1), 2);    
%                         
                          
%   GD_Handle=@(x,varargin)[cos(pi*x(:,2)) cos(pi*x(:,1))];
%     GD_Handle=@g_;
    
    QuadRule=Duffy(TProd(gauleg(0,1,10)));
    NREFS =;
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1 0;1 0;1 1;-1 1];
    Mesh.Elements = [1 2 4;2 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = [-1 -1 -2 -2];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    
    
    % Initialize mesh
%    DHANDLE = @dist_circ;                                    % Signed distance function
%    CEN = [0 0];                                               % Center of the circle
%    RAD = 1;                                                   % Radius of the circle
%    Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
%    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
%    Mesh = add_Edges(Mesh);
%    Mesh = add_Edge2Elem(Mesh);
%    Loc = get_BdEdges(Mesh);
%    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
%    Mesh.BdFlags(Loc) = -1;
   
   err=zeros(NREFS,1);
   h=zeros(NREFS,1);
   Dofs=zeros(NREFS,1);
   approx=zeros(NREFS,3);
    
    for i = 1:NREFS
       Mesh=refine_REG(Mesh);
%         Mesh = refine_REG(Mesh,DHANDLE,CEN,RAD);    
         Mesh = add_Edge2Elem(Mesh);  
         Mesh =add_VBdFlags(Mesh);
%         b_Mesh=refine_BAR(Mesh);
%         b_Mesh = add_Edge2Elem(b_Mesh);
%         b_Mesh=add_AdjNodes(b_Mesh);
%         TD=assemMat_TransOneD(Mesh,b_Mesh);
%         TP=assemMat_TransOneP(Mesh,b_Mesh);
        
        Mesh=add_Patches(Mesh);
        M = assemMat_W1F(Mesh,@MASS_W1F,MU_HANDLE, P7O6());
        sigma=SIGMA_Handle(Mesh.Coordinates);
        Mdis = assemMat_W1Fdis(Mesh,sigma);
        MassZero=assemMat_MassZeroD(Mesh);
        MassOne=assemMat_MassOneD(Mesh);
        MassTwo=assemMat_MassTwoD(Mesh);
        TopRot=assemMat_TopRot(Mesh);
        cU=assemCochain_1f(Mesh,U_EX_Handle,gauleg(0,1,10));
        cV=assemCochain_1f(Mesh,V_Handle,gauleg(0,1,10));
        
        cUV=assemCochainD_2f(Mesh,V_U_EX_Handle,QuadRule);
        cUVp=assemCochain_2f(Mesh,V_U_EX_Handle,QuadRule);
        cB=TopRot*cU;
        cUV2=assemCochain_0f(Mesh,V_U_EX_Handle);
       
        cUV2=MassZero*cUV2;
        sigma=SIGMA_Handle(Mesh.Coordinates);
        En=assemWedge_1f1f(Mesh,cU,cV,sigma);
        UxB=assemWedge_2f1f(Mesh,cB,cU);
%         plot_P0Cell(En,Mesh);
%         colorbar;
         plot_P0Cell(cUV,Mesh);
         colorbar;
%         plot_P0Cell(cUV2,Mesh);
%         colorbar;
         plot_P0Simplex(cUVp,Mesh);
         colorbar;
  %      L2Err_D2f(Mesh,En,QuadRule,V_U_EX_Handle)
        L2Err_D2f(Mesh,cUV,QuadRule,V_U_EX_Handle)
        L2Err_P2f(Mesh,cUVp,QuadRule,V_U_EX_Handle)
   end             
     
%          M_P0=assemMat_P0(Mesh,@MASS_P0,MU_HANDLE,P7O6());
%         Curl_U_L=assemLoad_P0(Mesh,P7O6(),CURL_U_EX_Handle);
%         curl_u=M_P0\Curl_U_L;
        
%         plot_PC(curl_u,Mesh);
%           colorbar;
%         
 
         L1=assemLoad_W1F(Mesh,P7O6(),SOL1_Handle);
         L2=assemLoad_W1F(Mesh,P7O6(),SOL2_Handle);
         L3=assemLoad_W1F(Mesh,P7O6(),SOL3_Handle);
         L4=assemLoad_W1F(Mesh,P7O6(),U_EX_Handle);
         L5=assemLoad_W1F(Mesh,P7O6(),Zero_Handle);
         
         L1RT=assemLoad_RT(Mesh,P7O6(),SOL1_Handle);
         L2RT=assemLoad_RT(Mesh,P7O6(),SOL2_Handle);
         
         ContrOne=assemMat_ContrOne(Mesh,V_Handle);
         ContrTwo=assemMat_ContrTwo(Mesh,V_Handle);
         ContrTwo2=assemMat_ContrTwo2(Mesh,V_Handle);

        % Mesh=add_UPWINDDATA(Mesh);
         
        % ContrRotUp=assemMat_W1F(Mesh,@CONTRROTUP, area, V_Handle, P706());
        ContrRot=assemMat_W1F(Mesh,@CONTRROT,V_Handle, P3O2());
        GradContr=assemMat_W1F(Mesh,@GRADCONTR,V_Handle, Grad_V_Handle,P3O2());
        TransOne=assemMat_W1F(Mesh,@TRANSONE,V_Handle, P3O2());
        
        TopGrad=assemMat_TopGrad(Mesh);
        TopRot=assemMat_TopRot(Mesh);
     
        MassOne=assemMat_MassOneD(Mesh); 
        MassTwo=assemMat_MassTwoD(Mesh);
%       MassZero=assemMat_MassZeroD(Mesh);
  
        A=M*ContrTwo2*TopRot;
        B=M*TopGrad*ContrOne;
        C=TopRot'*MassTwo*TopRot;
     
       R=assemMat_W1F(Mesh,@STIMA_Curl_W1F,EPSI_Handle,P3O3());
         
%          plot_Norm_W1F(TopRot'*MassTwo*TopRot*U_ex,Mesh);
%          colorbar;
        plot_Norm_W1F(U_ex,Mesh);
        colorbar;
%           figure;
%           plot_W1F(M*ContrTwo*TopRot*U_ex,Mesh);
%           figure;
      %     plot_Norm_W1F(GradContr*U_ex,Mesh);
%            colorbar;
%            plot_Norm_W1F(M\L2,Mesh);
%            colorbar;         
%           plot_Norm_W1F(GradContr*U_ex,Mesh);
%           colorbar;
%          plot_Norm_W1F(M*ContrTwo*TopRot*U_ex,Mesh);
%          colorbar;
%          figure;
%           plot_LFE(ContrOne*U_ex,Mesh);
%           colorbar;
%         % System Matrix
     
   % Direchlet Boundary data
     [U,FreeDofs] = assemDir_W1F(Mesh,[-1],U_EX_Handle,gauleg(0,1,10));

   % Load vector
      L=L2+L3;
      S=B+C;
      
      L = L - S*U;
      % Solve the system
     
      U(FreeDofs) = S(FreeDofs,FreeDofs)\L(FreeDofs);
      err(i)=L2Err_W1F(Mesh,U,QuadRule,U_EX_Handle)
      h(i)=get_MeshWidth(Mesh);
      Dofs(i)=size(Mesh.Edges,1);
      plot_Norm_W1F(U,Mesh);
      colorbar; 
    %end
    
    fig = figure('Name','Discretization error');
    plot(Dofs,err,'ro--'); grid('on');
    set(gca,'XScale','log','YScale','log');
    xlabel('{\bf Dofs}');
    ylabel('{\bf Error}');
    p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(NREFS-3:NREFS)),1);
    add_Slope(gca,'East',p(1),'r-');
    % Clear memory'
    
 clear all
    