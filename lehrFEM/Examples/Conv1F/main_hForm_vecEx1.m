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
JIG=1;                   % 
NREFS =5;             % number of refinements
QuadRule = P7O6();

% model parameters
isigma = 1;             % diffusion constant
mu = 1                  %
supg1=1;                 %impact of supg modification
supg2=1;                 %impact of supg modification

% select test code
d=getData_EsHvEx1();

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
    ContrOne=assemMat_Contr1f(NewMesh,@(x)-d.V_Handle(x));  % contraction of one forms
    V=d.V_Handle(NewMesh.Coordinates);
    ContrTwo=assemMat_Contr2f(NewMesh,-V);   % contraction of two forms
    
    % topological derivatives
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);      % topological Rotation

    % (-v x curl u, u'): standard local quadrature
    ContrRot = assemMat_W1F(NewMesh,@STIMA_ContrRot,d.V_Handle, QuadRule);            % -v x rot u FEM
    
    LieInnM = -assemMat_Lie_Inn_W1F(NewMesh,@(x,t)-d.V_Handle(x));
    LieInnP = assemMat_Lie_Inn_W1F(NewMesh,@(x,t)d.V_Handle(x));
    LieInn=1/2*(LieInnM+LieInnP);
    LieVol= assemMat_W1F(NewMesh,@STIMA_ContrRot,d.V_Handle, P7O6());
    LieVol = LieVol+ assemMat_W1F(NewMesh,@STIMA_GradContr,d.V_Handle, gauleg(0,1,3));
  
    % upwind
    DD = TopRot'*MassTwo*TopRot;             % curl curl u geom.
    ID_up = M_W1F*ContrTwo*TopRot;          % -v x curl u geom.
    DI_up =  M_W1F*TopGrad*ContrOne;      % grad(v.u) geom.
    
    %  stiffnes matrices
    A_st = (isigma*DD-(LieInn'+LieVol')+M_W1F);
    A_up = (isigma*DD+ID_up'+DI_up'+M_W1F);
    % A_supg =(A_st+D_supg);

    % source term
    L = assemLoad_W1F(NewMesh,P1O2(),d.CurlFs_Handle,isigma); 
    %L_supg=assemLoad_W1F_SUPG(NewMesh,P1O2(),d.V_Handle,d.CurlFs_Handle,isigma,supg1,supg2);
    
    % Direchlet boundary
    [Hv_up,FreeDofs] = assemDir_W1F(NewMesh,-1,d.Hv_Handle,gauleg(0,1,10),isigma);
    Hv_st = Hv_up;
    
    % Rot Contr Boundary Contribution
    L_up = L- A_up*Hv_up;
    L_st = L - A_st*Hv_st;
    
    % solving system
    Hv_up(FreeDofs) = A_up(FreeDofs,FreeDofs)\L_up(FreeDofs);
    Hv_st(FreeDofs) = A_st(FreeDofs,FreeDofs)\L_st(FreeDofs);
    
    % calulating L^2-error
    L2err(i,1) = L2Err_W1F(NewMesh,Hv_up,P7O6(),d.Hv_Handle,0,isigma);
    L2err(i,2) = L2Err_W1F(NewMesh,Hv_st,P7O6(),d.Hv_Handle,0,isigma); 
    
    % calculating H@1-semi -error
    HCurlSerr(i,1) = HCurlSErr_W1F(NewMesh,Hv_up,P7O6(),d.CURL_Hv_Handle,0,isigma);
    HCurlSerr(i,2) = HCurlSErr_W1F(NewMesh,Hv_st,P7O6(),d.CURL_Hv_Handle,0,isigma);
    
    Dofs(i)=size(NewMesh.Coordinates,1);
     
%    plot_Norm_W1F(Hv_up,NewMesh); colorbar; title('Upwind');
%     plot_Norm_W1F(Hv_st,NewMesh); colorbar; title('Standard');

     plot_P0(TopRot*Hv_up,NewMesh); colorbar; title('Upwind');
     plot_P0(TopRot*Hv_st,NewMesh); colorbar; title('Standard');
    
end
fig = figure('Name','L^2error Hv');
plot(h,h*(L2err(1,1)/h(1)),'y--',h,L2err(:,1),'r-',h,L2err(:,2),'b-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf L^2-Error}');
% legend('h','Upwind','Standard','UpwindC','UpwindVertex','StandardUp','Location','NorthEast')
legend('h','Upwind','Standard')
%p = polyfit(log(h(1:NREFS)),log(L2err(1:NREFS,1)),1);
%add_Slope(gca,'SouthEast',p(1),'r--');
%p = polyfit(log(h(1:NREFS)),log(L2err(1:NREFS,2)),1);
%add_Slope(gca,'East',p(1),'b-.');
%p = polyfit(log(h(1:NREFS)),log(L2err(1:NREFS,3)),1);
%add_Slope(gca,'South',p(1),'g-');

fig = figure('Name','HCurl-semi-error v');
plot(h,h*(HCurlSerr(1,1)/h(1)),'y--',h,HCurlSerr(:,1),'r-',h,HCurlSerr(:,2),'b-'); grid('on');
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
