% Run script discrete differential forms.
% h formulation, scalar component

%   Copyright 2007-2009 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


% driver routine for several diffusion convection discretizations
% see below

clear Mesh
% Initialize constants
JIG=2;                   % 
NREFS =6;             % number of refinements

% model parameters
isigma = 10^-9;             % diffusion constant
mu = 1                  %
supg1=1;                 %impact of supg modification
supg2=1;                 %impact of supg modification

% select test code
d=getData_EvHs();

QuadRule = P7O6();

Mesh.Coordinates =d.Coordinates;
Mesh.Elements = d.Elements;
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = d.boundtype;
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);

% allocate memory
L2err=zeros(NREFS,4);
H1Serr=zeros(NREFS,4);
h=zeros(NREFS,1);
Dofs=zeros(NREFS,1);

for i = 1:NREFS

    %refine Mesh
    Mesh = refine_REG(Mesh);

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
    
    h(i)=get_MeshWidth(NewMesh);

    %Laplace
    A = assemMat_LFE(NewMesh,@STIMA_Lapl_LFE);

    %Mass
    M = assemMat_LFE(NewMesh,@MASS_LFE);
    MassZero=assemMat_Mass0fD(NewMesh);          %diagonal Mass Matrix
    
    %Convection term 
       
    % upwinding 
    %           D_up=assemMat_LFE(NewMesh, @STIMA_ContrGrad_Up, d.V_Handle);
    %           or equivalently
    TopRot=assemMat_TopRot(NewMesh);      % topological Rotation
    TopGrad=assemMat_TopGrad(NewMesh);
    ContrOne=assemMat_Contr1f(NewMesh,@(x)-d.V_Handle(x));
    D_up=MassZero*ContrOne*TopGrad;
    %            D_up=M*ContrOne*TopGrad;

    % standard Galerkin
    D_st=assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,d.V_Handle,P1O2());

    % SUPG - stabelization
    D_supg = assemMat_LFE(NewMesh,@STIMA_SUPG_LFE,P1O2(), d.V_Handle,isigma,supg1,supg2);
  
    %  stiffnes matrices
    A_up =(isigma*A-D_up'+M);
    A_st =(isigma*A+D_st'+M);
    A_supg =(A_st-D_supg);

    % source term
    L = assemLoad_LFE(NewMesh,P1O2(),d.CurlFv_Handle,isigma);   
    L_supg=assemLoad_LFE_SUPG(NewMesh,P1O2(),d.V_Handle,d.CurlFv_Handle,isigma,supg1,supg2);
    
    % Direchlet boundary
    [Hs_up,FreeDofs] = assemDir_LFE(NewMesh,-1,d.Hs_Handle,isigma);
    Hs_st=Hs_up;
    Hs_supg=Hs_up;
  
    L_up = -L - A_up*Hs_up;
    L_st = -L - A_st*Hs_st;
    L_supg = -(L+L_supg) - A_supg*Hs_supg;

    % solving system
    Hs_up(FreeDofs) = A_up(FreeDofs,FreeDofs)\L_up(FreeDofs);
    Hs_st(FreeDofs) = A_st(FreeDofs,FreeDofs)\L_st(FreeDofs);
    Hs_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);

    % calulating L^2-error
    L2err(i,1) = L2Err_LFE(NewMesh,Hs_up,P7O6(),d.Hs_Handle,0,isigma);
    L2err(i,2) = L2Err_LFE(NewMesh,Hs_st,P7O6(),d.Hs_Handle,0,isigma);
    L2err(i,3) = L2Err_LFE(NewMesh,Hs_supg,P7O6(),d.Hs_Handle,0,isigma)
    
    % calculating H@1-semi -error
    H1Serr(i,1) = H1SErr_LFE(NewMesh,Hs_up,P7O6(),d.GRAD_Hs_Handle,0,isigma);
    H1Serr(i,2) = H1SErr_LFE(NewMesh,Hs_st,P7O6(),d.GRAD_Hs_Handle,0,isigma);
    H1Serr(i,3) = H1SErr_LFE(NewMesh,Hs_supg,P7O6(),d.GRAD_Hs_Handle,0,isigma);

    Dofs(i)=size(NewMesh.Coordinates,1);
     
    plot_LFE(Hs_up,NewMesh); colorbar; title('Upwind');
    plot_LFE(Hs_st,NewMesh);  colorbar; title('Standard');
    plot_LFE(Hs_supg,NewMesh); colorbar;('SUPG');
    %plot_Mesh(NewMesh,'a'); title('mesh');
    %L_H = assemLoad_LFE(NewMesh,P3O3(),d.Es_Handle,isigma);
    %H=M\L_H;
    %plot_LFE(H,NewMesh);colorbar;
    
    %plot_Norm_W1F(TopGrad*Es_up,NewMesh); colorbar; title('Upwind');
    %plot_Norm_W1F(TopGrad*Es_st,NewMesh);  colorbar; title('Standard');
    %plot_Norm_W1F(TopGrad*Es_supg,NewMesh); colorbar;('SUPG');
        
end;

% exact solution
% L_H = assemLoad_LFE(NewMesh,P1O2(),d.H_EX_Handle,a);
% H(FreeDofs)=M(FreeDofs,FreeDsofs)\L_H(FreeDofs);
% plot_LFE(H,NewMesh);
% plot_Mesh(NewMesh,'a');

% plot L2 error;
fig = figure('Name','L^2error Hs');
plot(h,h*(L2err(1,1)/h(1)),'c--',h,L2err(:,1),'r-',h,L2err(:,2),'b-',h,L2err(:,3),'g-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf L^2-Error}');
legend('h','Upwind','Standard','SUPG','Location','NorthEast')
p = polyfit(log(h(1:NREFS)),log(L2err(1:NREFS,1)),1);
add_Slope(gca,'SouthEast',p(1),'r-');
p = polyfit(log(h(1:NREFS)),log(L2err(1:NREFS,2)),1);
add_Slope(gca,'East',p(1),'b-');
p = polyfit(log(h(1:NREFS)),log(L2err(1:NREFS,3)),1);
add_Slope(gca,'South',p(1),'g-');

fig = figure('Name','H^1-semi-error Hs');
plot(h,h*(H1Serr(1,1)/h(1)),'c--',h,H1Serr(:,1),'r-',h,H1Serr(:,2),'b-',h,H1Serr(:,3),'g-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf L^2-Error}');
legend('h','Upwind','Standard','SUPG','Location','NorthEast')
p = polyfit(log(h(1:NREFS)),log(H1Serr(1:NREFS,1)),1);
add_Slope(gca,'SouthEast',p(1),'r-');
p = polyfit(log(h(1:NREFS)),log(H1Serr(1:NREFS,2)),1);
add_Slope(gca,'East',p(1),'b-');
p = polyfit(log(h(1:NREFS)),log(H1Serr(1:NREFS,3)),1);
add_Slope(gca,'South',p(1),'g-');