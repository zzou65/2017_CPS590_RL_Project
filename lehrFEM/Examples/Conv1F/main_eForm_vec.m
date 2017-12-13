% Run script discrete differential forms.
% main h formulation vetorial component

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

clear Mesh;
% close all;
% clear all;

% Initialize constants
JIG=1;                        % 
NREFS =7;                  % number of refinements
QuadRule = P7O6();

% model parameters
isigma = 10^-8;            % diffusion constant
mu = 1                        %
supg1=1;                     %impact of supg modification
supg2=1;                     %impact of supg modification

% select test code
d=getData_EvHs();

% Initialize mesh
Mesh.Coordinates = d.Coordinates;
Mesh.Elements = d.Elements;
Mesh = add_Edges(Mesh);

Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = d.boundtype;
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);

L2err=zeros(NREFS,4);
HCurlSerr=zeros(NREFS,4);
h=zeros(NREFS,1);
cond(NREFS,2);
Dofs=zeros(NREFS,1);

for i = 1:NREFS
    i

    Mesh=refine_REG(Mesh);
    Mesh=add_Edge2Elem(Mesh);

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

    % Meshwidth and DOFS
    h(i)=get_MeshWidth(NewMesh);
    Dofs(i)=size(NewMesh.Edges,1);
    
    % several mass matrices
    % 0-forms
    %M_LFE=assemMat_LFE(NewMesh,@MASS_LFE);
    %MassZero=assemMat_Mass0fD(NewMesh);         %diagonal 
    % 1-forms
    M_W1F = assemMat_W1F(NewMesh,@MASS_W1F,@(x,vargin)1, P7O6());            % Mass Matrix Whitney-1-Forms
    %MassOne=assemMat_Mass1fD(NewMesh);                                                  % diagonal Mass Matrix of 1-Forms
    % 2-forms
    MassTwo=assemMat_Mass2fD(NewMesh);                                                    % diagonal Mass Matrix of Two-Forms

    % contraction
    ContrOne=assemMat_Contr1f(NewMesh,@(x)d.V_Handle(x));  % contraction of one forms
    V=d.V_Handle(NewMesh.Coordinates);
    ContrTwo=assemMat_Contr2f(NewMesh,V);   % contraction of two forms
    ContrTwoV=assemMat_Contr_1FVertUP(NewMesh,@(x)d.V_Handle(x));   % contraction of two forms( vertexwise upwind)
    
    % topological derivatives
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);      % topological Rotation

    % SUPG-stabelization  ?????
    %D_supg = assemMat_W1F(NewMesh,@STIMA_SUPG_W1F, d.V_Handle, QuadRule);

    % (-v x curl u, u'): standard local quadrature
    ContrRot = assemMat_W1F(NewMesh,@STIMA_ContrRot,d.V_Handle, QuadRule);            % -v x rot u FEM
    ContrRotUp = assemMat_ContrRot_UPQuad(NewMesh,@(x)d.V_Handle(x));
    % (grad( v u).u'): mixed ansatz (grad(k),u') with (k,q)=(v u q), q, k 0-forms 
    %ContrOne_fem=assemMat_Contr1f_FEM(NewMesh,d.V_Handle,P7O6());    %contraction of one forms, exact quadrature
    %GradContr=M_W1F*TopGrad*(M_LFE\ContrOne_fem);
    %ContrOne_fem=assemMat_Contr1f_FEM(NewMesh,d.V_Handle,P3O2());      % contraction of one forms, non-exact quadrature 
    %GradContr=M_W1F*TopGrad*(MassZero\ContrOne_fem);
    
    %LieInn_FEM = assemMat_Lie_Inn_W1F(NewMesh,d.V_Handle);
    %LieVol_FEM = assemMat_W1F(New_Mesh,@STIMA_Lie,d.V_Handle,d.JacV_Handle);
  
    % Stiffness matrices, 
    % upwind
    DD = TopRot'*MassTwo*TopRot;             % curl curl u geom.
    ID_up = ContrTwo*TopRot;                     % -v x curl u geom.
    ID_upV = ContrTwoV*TopRot;
    ID_st = ContrRot;
    ID_stUp = ContrRotUp;
    
    DI_up = TopGrad*ContrOne;                   % grad(v.u) geom.
    %ID_st = GradContr;                               % grad(v.u) geom.
    %ID_up = GradContr;                              % grad(v.u) geom.
    
    %  stiffnes matrices
    A_up =(isigma*DD+M_W1F*ID_up+M_W1F);
    A_upV =(isigma*DD+M_W1F*ID_upV+M_W1F);
    A_st =(isigma*DD+ID_st+M_W1F);
    A_stUp =(isigma*DD+ID_stUp+M_W1F);
    A_upc =(isigma*DD+M_W1F*ID_up+M_W1F*DI_up+M_W1F);
    % A_supg =(A_st+D_supg);

    % source term
    L = assemLoad_W1F(NewMesh,P7O6(),d.Fv_Handle,isigma); 
    %L_supg=assemLoad_W1F_SUPG(NewMesh,P1O2(),d.V_Handle,d.CurlFs_Handle,isigma,supg1,supg2);
    
    % Direchlet boundary
    [Ev_up,FreeDofs] = assemDir_W1F(NewMesh,[-1 -2],d.Ev_Handle,gauleg(0,1,10),isigma);
    Ev_st = Ev_up;
    Ev_upV = Ev_up;
    Ev_stUp = Ev_up;
    Ev_upc = Ev_up;
    %Hv_supg = Hv_up;
    
    % Rot Contr Boundary Contribution
    %L_bnd = assemRotContr_Bnd_W1F(NewMesh, -1, d.Ev_Handle, d.V_Handle, gauleg(0,1,10),isigma);
    c1Bd=assemBndCochain_1f(NewMesh,[-1 -2],d.VCURL_Ev_Handle,gauleg(0,1,10),isigma);
    %cVUBd=assemBndCochain_0f(NewMesh,[-2 -3],d.V_Ev_Handle,isigma);  % boundary data on inflow boundary

    L_up = L - A_up*Ev_up;%-M_W1F*c1Bd;
    L_upV = L - A_upV*Ev_upV;%--M_W1F*c1Bd;
    L_st = L - A_st*Ev_st;
    L_stUp = L - A_stUp*Ev_stUp;%--M_W1F*c1Bd;
    L_upc = L - A_upc*Ev_upc;% + M_W1F*TopGrad*cVUBd+M_W1F*c1Bd;%;
    %L_supg = -(L+L_supg) - A_supg*Hv_supg;
    
    % solving system
    Ev_up(FreeDofs) = A_up(FreeDofs,FreeDofs)\L_up(FreeDofs);
    Ev_upV(FreeDofs) = A_upV(FreeDofs,FreeDofs)\L_upV(FreeDofs);
    Ev_st(FreeDofs) = A_st(FreeDofs,FreeDofs)\L_st(FreeDofs);
    Ev_stUp(FreeDofs) = A_stUp(FreeDofs,FreeDofs)\L_stUp(FreeDofs);
    Ev_upc(FreeDofs) = A_upc(FreeDofs,FreeDofs)\L_upc(FreeDofs);
    
    cond(i,1) = condest(A_up(FreeDofs,FreeDofs));
    cond(i,2) = condest(A_up(FreeDofs,FreeDofs));
    
    %Hv_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);
    %Ev_upc=Ev_upc-DI_up*Ev_upc;
     
    % calulating L^2-error
    L2err(i,1) = L2Err_W1F(NewMesh,Ev_up,P7O6(),d.Ev_Handle,0,isigma);
    L2err(i,2) = L2Err_W1F(NewMesh,Ev_st,P7O6(),d.Ev_Handle,0,isigma); 
    L2err(i,3) = L2Err_W1F(NewMesh,Ev_upc-DI_up*Ev_upc,P7O6(),d.Ev_Handle,0,isigma);
    L2err(i,4) = L2Err_W1F(NewMesh,Ev_upV,P7O6(),d.Ev_Handle,0,isigma);
    L2err(i,5) = L2Err_W1F(NewMesh,Ev_stUp,P7O6(),d.Ev_Handle,0,isigma)
    %L2err(i,3) = L2Err_W1F(NewMesh,Hv_supg,P7O6(),d.Hv_Handle,0,isigma)
    
    % calculating H@1-semi -error
    HCurlSerr(i,1) = HCurlSErr_W1F(NewMesh,Ev_up,P7O6(),d.CURL_Ev_Handle,0,isigma);
    HCurlSerr(i,2) = HCurlSErr_W1F(NewMesh,Ev_st,P7O6(),d.CURL_Ev_Handle,0,isigma);
    HCurlSerr(i,3) = HCurlSErr_W1F(NewMesh,Ev_upc+DI_up*Ev_upc,P7O6(),d.CURL_Ev_Handle,0,isigma);
    HCurlSerr(i,4) = HCurlSErr_W1F(NewMesh,Ev_upV,P7O6(),d.CURL_Ev_Handle,0,isigma);
    HCurlSerr(i,5) = HCurlSErr_W1F(NewMesh,Ev_stUp,P7O6(),d.CURL_Ev_Handle,0,isigma);
   
    %HCurlSerr(i,3) = HCurlSErr_W1F(NewMesh,Hv_supg,P7O6(),d.CURL_Hv_Handle,0,isigma);

    Dofs(i)=size(NewMesh.Coordinates,1);
     
     plot_Norm_W1F(Ev_upV,NewMesh); colorbar; title('Upwind');
     plot_Norm_W1F(Ev_st,NewMesh); colorbar; title('Standard');
%     plot_Norm_W1F(Ev_upV,NewMesh); colorbar; title('UpwindVertex');
%     plot_Norm_W1F(Ev_stUp,NewMesh); colorbar; title('StandardUp');

%     plot_P0(TopRot*Ev_up,NewMesh); colorbar; title('Upwind');
%     plot_P0(TopRot*Ev_st,NewMesh); colorbar; title('Standard');
%     plot_P0(TopRot*Ev_upV,NewMesh); colorbar; title('UpwindVertex');
%     plot_P0(TopRot*Ev_stUp,NewMesh); colorbar; title('StandardUp');
    
    %plot_Norm_W1F(Hv_supg,NewMesh); colorbar;('SUPG');
    %figure; plot_W1F(Ev_up,NewMesh); title('Upwind');
    %figure; plot_W1F(Ev_upc,NewMesh);  title('Standard');

    %plot_Mesh(NewMesh,'a'); title('mesh');
%     L_H = assemLoad_W1F(NewMesh,P3O3(),d.Ev_Handle,isigma);
%     H = M_W1F\L_H;
    
    %  plot_Norm_W1F(A_upc*H,NewMesh);colorbar;
    %  plot_Norm_W1F(L,NewMesh); colorbar;
    %  plot_Norm_W1F(A_up*H-DI_up*H,NewMesh);colorbar;
    %  plot_Norm_W1F(L,NewMesh);colorbar;

    %plot_Norm_W1F(H,NewMesh);colorbar; 
    %plot_P0(TopRot*H,NewMesh); colorbar; title('Analytic');
    %plot_P0(TopRot*Ev_up,NewMesh); colorbar; title('Upwind');
    %plot_P0(TopRot*Ev_upc,NewMesh); colorbar; title('UpwindC');
    %[norm(Ev_upc-Ev_up-DI_up*Ev_up) norm(ContrTwo*TopRot*(Ev_up-Ev_upc))]
    
    %[condest(A_up(FreeDofs,FreeDofs)) condest(A_upc(FreeDofs,FreeDofs))]

        
end
fig = figure('Name','L^2error Ev');
plot(h,h.^(1/2)*(L2err(1,3)/h(1)),'y--',h,L2err(:,1),'r-',h,L2err(:,2),'b-',h,L2err(:,3),'g-',h,L2err(:,4),'k-',h,L2err(:,5),'c-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf L^2-Error}');
legend('h','Upwind','Standard','UpwindC','UpwindVertex','StandardUp','Location','NorthEast')
%p = polyfit(log(h(1:NREFS)),log(L2err(1:NREFS,1)),1);
%add_Slope(gca,'SouthEast',p(1),'r--');
%p = polyfit(log(h(1:NREFS)),log(L2err(1:NREFS,2)),1);
%add_Slope(gca,'East',p(1),'b-.');
%p = polyfit(log(h(1:NREFS)),log(L2err(1:NREFS,3)),1);
%add_Slope(gca,'South',p(1),'g-');

fig = figure('Name','HCurl-semi-error Ev');
plot(h,h*(HCurlSerr(1,1)/h(1)),'y--',h,HCurlSerr(:,1),'r-',h,HCurlSerr(:,2),'b-',h,HCurlSerr(:,3),'g-',h,HCurlSerr(:,4),'k-',h,HCurlSerr(:,5),'c-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf LCurl-semi-Error}');
legend('h','Upwind','Standard','UpwindC','UpwindVertex','StandardUp','Location','NorthEast')
% p = polyfit(log(h(1:NREFS)),log(HCurlSerr(1:NREFS,1)),1);
% add_Slope(gca,'SouthEast',p(1),'r--');
% p = polyfit(log(h(1:NREFS)),log(HCurlSerr(1:NREFS,2)),1);
% add_Slope(gca,'East',p(1),'b-.');
% p = polyfit(log(h(1:NREFS)),log(HCurlSerr(1:NREFS,3)),1);
% add_Slope(gca,'South',p(1),'g-');
