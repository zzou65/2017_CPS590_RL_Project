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
e_a=1;                           % -(v x curl u)
e_b=1;  e_bf=1;            % (grad v.u)
e_c=10^-8;                   % curlcurl u
e_m =1;                         % u
d=getData(15);             % struct containing Coordinates, Elements and various function handles
EPSI_Handle = @(x,varargin)ones(size(x,1),1);
MU_HANDLE = @(x,varargin)1;
QuadRule=Duffy(TProd(gauleg(0,1,2)));
NREFS = 6;
JIG =2;

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
err5=zeros(NREFS,1);
err6=zeros(NREFS,1);
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
    New_Mesh=add_Edge2Elem(New_Mesh);
    New_Mesh = add_DGData(New_Mesh);
    % contraction
    ContrOne=assemMat_Contr1f(New_Mesh,d.V_Handle);  % contraction of one forms

    V=d.V_Handle(New_Mesh.Coordinates);
    ContrTwo=assemMat_Contr2f(New_Mesh,V);   % contraction of two forms

    % topological derivatives
    TopGrad=assemMat_TopGrad(New_Mesh);   % topological Gradient
    TopRot=assemMat_TopRot(New_Mesh);      % topological Rotation

    % several mass matrices
    % 0-forms
    M_LFE=assemMat_LFE(New_Mesh,@MASS_LFE);
    MassZero=assemMat_Mass0fD(New_Mesh);         %diagonal 
    % 1-forms
    M_W1F = assemMat_W1F(New_Mesh,@MASS_W1F,MU_HANDLE, P7O6());  % Mass Matrix Whitney-1-Forms
    MassOne=assemMat_Mass1fD(New_Mesh);                                          % diagonal Mass Matrix of 1-Forms
    % 2-forms
    MassTwo=assemMat_Mass2fD(New_Mesh);                                          % diagonal Mass Matrix of Two-Forms

    % Stiffness matrices,
    % standard
    %LieInnM = -assemMat_Lie_Inn_W1F(New_Mesh,@(x,t)-d.V_Handle(x));
    %LieInnP = assemMat_Lie_Inn_W1F(New_Mesh,@(x,t)d.V_Handle(x));
    %LieInn_FEM=1/2*(LieInnM+LieInnP);
    %LieVol_FEM= assemMat_W1F(New_Mesh,@STIMA_ContrRot,d.V_Handle, P7O6());
    %LieVol_FEM= LieVol_FEM+ assemMat_W1F(New_Mesh,@STIMA_GradContr,d.V_Handle, gauleg(0,1,3));
    Lie1 = assemMat_W1F(New_Mesh,@STIMA_ContrRot,d.V_Handle, P7O6());
    Lie2 =  assemMat_Inn_LieW1F(New_Mesh,@STIMA_Inn_LieW1F,d.V_Handle,QuadRule2D);

    % Stiffness matrices, 
    % upwind
    A=M_W1F*ContrTwo*TopRot;               % -v x curl u geom.
    B=M_W1F*TopGrad*ContrOne;             % grad(v.u) geom.
    C=TopRot'*MassTwo*TopRot;               % curl curl u geom.
  
    % boundary data on inflow boundary and righthandside
    c1=assemCochain_1f(New_Mesh,d.SOL1_Handle,gauleg(0,1,10));     % -v x curl u righthand side
    c2=assemCochain_1f(New_Mesh,d.SOL2_Handle,gauleg(0,1,10));     % grad(v.u) righthand side
    c3=assemCochain_1f(New_Mesh,d.SOL3_Handle,gauleg(0,1,10));     % curl curl rightand side
    c4=assemCochain_1f(New_Mesh,d.U_EX_Handle,gauleg(0,1,10));      % u righthand side
    cVUBd=assemBndCochain_0f(New_Mesh,[-2 -3],d.V_U_EX_Handle);  % boundary data on inflow boundary
    c1Bd=assemBndCochain_1f(New_Mesh,[-2],d.SOL1_Handle,gauleg(0,1,10)); % boundary

    % exact solution
    cU=assemCochain_1f(New_Mesh,d.U_EX_Handle,gauleg(0,1,10));

    % calculate solutions
    % free DOFS
    nonInFlow=find(New_Mesh.BdFlags ~= -1);
    nonOutFlow=find(New_Mesh.BdFlags ~= -2);
    FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
    
    % upwind scheme
    S_up=e_m*M_W1F+e_a*A+e_c*C+e_b*B;
    cf=e_m*c4+e_a*c1+e_c*c3+e_bf*c2;
    %h1=c1;
    %h2=c2;
    %FreeDofs=nonInFlow;
    %h1(nonInFlow)=0;
    %FreeDofs=nonOutFlow;
    cU_up=cU;
    cU_up(FreeDofs)=0;
    %Lf_up=M_W1F*(cf-e_a*h1)-S_up*cU_up;
    Lf_up=M_W1F*(cf)-S_up*cU_up;

    % solve system
    cU_up(FreeDofs)=(S_up(FreeDofs,FreeDofs))\Lf_up(FreeDofs);
    norm(S_up(FreeDofs,FreeDofs)*cU_up(FreeDofs)-Lf_up(FreeDofs))

    % FEM-based upwind scheme
    S_upFEM=e_m*M_W1F+e_c*C+e_a*Lie1+e_b*Lie2;
    cf=e_m*c4+e_a*c1+e_c*c3+e_bf*c2;
    %h1=c1;
    %h2=c2;
    %FreeDofs=nonInFlow;
    %h1(nonInFlow)=0;F
    %FreeDofs=nonOutFlow;
    cU_upFEM=cU;
    cU_upFEM(FreeDofs)=0;
    %Lf_up=M_W1F*(cf-e_a*h1)-S_up*cU_up;
    Lf_upFEM=M_W1F*cf-S_upFEM*cU_upFEM;

    % solve system
    cU_upFEM(FreeDofs)=S_upFEM(FreeDofs,FreeDofs)\Lf_upFEM(FreeDofs);
 
    %subplot(2,2,3); 
    %plot_Norm_W1F(cU_up,New_Mesh); colorbar; title('U Interpol');
    %plot_P0(TopRot*cU_up,New_Mesh); colorbar; title('U Interpol');
    %subplot(2,2,4); 
    %plot_Norm_W1F(cU_upFEM,New_Mesh); colorbar; title('U Patches');
    %plot_P0(TopRot*cU_upFEM,New_Mesh); colorbar; title('U Patches');
     
    err3(i)=L2Err_W1F(New_Mesh,cU_up,QuadRule,d.U_EX_Handle);
    err4(i)=HCurlSErr_W1F(New_Mesh,cU_up,QuadRule,d.CURL_U_EX_Handle);
    err5(i)=L2Err_W1F(New_Mesh,cU_upFEM,QuadRule,d.U_EX_Handle);
    err6(i)=HCurlSErr_W1F(New_Mesh,cU_upFEM,QuadRule,d.CURL_U_EX_Handle);

    h(i)=get_MeshWidth(New_Mesh);
    Dofs(i)=size(New_Mesh.Edges,1);

   % [min(1/2*eig(full(Ss(FreeDofs,FreeDofs)+Ss(FreeDofs,FreeDofs)'))) ...
   %  min(1/2*eig(full(S_up(FreeDofs,FreeDofs)+S_up(FreeDofs,FreeDofs)')))]
end

fig = figure('Name','Discretization error');
plot(Dofs,err3,'bo-',Dofs,err4,'b--',Dofs,err5,'co-',Dofs,err6,'c--'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf Dofs}');
ylabel('{\bf Error}');
p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err1(NREFS-3:NREFS)),1);
add_Slope(gca,'SouthWest',p(1),'r-');
p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS)),1);
add_Slope(gca,'SouthEast',p(1),'r--');
p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err3(NREFS-3:NREFS)),1);
add_Slope(gca,'SouthWest',p(1),'b-');
p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err4(NREFS-3:NREFS)),1);
add_Slope(gca,'SouthEast',p(1),'b--');
p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err5(NREFS-3:NREFS)),1);
add_Slope(gca,'SouthWest',p(1),'c-');
p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err6(NREFS-3:NREFS)),1);
add_Slope(gca,'SouthEast',p(1),'c--');


legend('L^2','Hcurl-semi','L^2(up)','Hcurl-semi(up)','L^2(upF)','Hcurl-semi(upF)','Location','NorthEast')
%clear all

fig = figure('Name','Discretization error');
plot(h,h*(err4(1)+err6(1))/(2*h(1)),'k--',h,err4,'b-',h,err6,'c-','linewidth', 3); grid('on');
set(gca,'YScale','log','XScale','log','FontWeight','bold');
xlabel('{\bf h}');
ylabel('{\bf Hcurl-semi error}');
legend('h','Upwind','Standard','Location','NorthEast')
%clear all