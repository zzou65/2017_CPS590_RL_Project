function main_standardSharp(sigma)
% Run script discrete differential forms.

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

clear Mesh;
filename= ['./results/exampleSharp_',num2str(sigma)];
tic;
% close all;
% clear all;

onlyStandard = 1;

% Initialize constant
% e_a* (v x curl u)+e_b*(grad v.u)+e_c* curlcurl u +e_m*u
e_a=1;                                            % -(v x curl u)
e_b=1;  e_bf=1; e_bfB=1;            % (grad v.u)
e_c=0;                                            % curlcurl u
e_m =0.1;                                         % u
d=getData(15);                               % struct containing Coordinates, Elements and various function handles
EPSI_Handle = @(x,varargin)ones(size(x,1),1);
MU_HANDLE=@(x,varargin)1;
QuadRule=Duffy(TProd(gauleg(0,1,10)));
NREFS =6;
JIG = 6;
%sigma=3/4;

% QuadRule2D.x=[0;1];
% QuadRule2D.w=[0.5;0.5];
QuadRule2D=gauleg(0,1,10);
% Initialize mesh

Mesh.Coordinates = d.Coordinates;
Mesh.Elements = d.Elements;
Mesh = add_Edges(Mesh);
Mesh = add_Edge2Elem(Mesh);

Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = d.boundtype;
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);

L2err1=zeros(NREFS,1);
L2err2=zeros(NREFS,1);
L2err3=zeros(NREFS,1);
Hcurlerr1=zeros(NREFS,1);
Hcurlerr2=zeros(NREFS,1);
Hcurlerr3=zeros(NREFS,1);
IDPen1=zeros(NREFS,1);
IDPen2=zeros(NREFS,1);
IDPen3=zeros(NREFS,1);
PIDPen1=zeros(NREFS,1);
PIDPen2=zeros(NREFS,1);
PIDPen3=zeros(NREFS,1);
JPen1=zeros(NREFS,1);
JPen2=zeros(NREFS,1);
JPen3=zeros(NREFS,1);
h=zeros(NREFS,1);
Dofs=zeros(NREFS,1);
approx=zeros(NREFS,3);
for i = 1:NREFS

    Mesh=refine_REG(Mesh);
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
        case 4
            nC=size(Mesh.Coordinates,1);
            Loc = [1:6,8:nC];
            FixedPos = zeros(nC,1);
            FixedPos(Loc) = 1;
            New_Mesh = jiggle(Mesh,FixedPos);
        case 5
            Loc = get_BdEdges(Mesh);
            Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
            Loc = nonzeros(Mesh.Vert2Edge(Loc,:));
            Loc = unique(Loc);
            Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
            FixedPos = zeros(size(Mesh.Coordinates,1),1);
            FixedPos(Loc) = 1;
            New_Mesh = jiggle(Mesh,FixedPos);
        case 6
            nSegments = 2.^floor((i)*sigma);
            points=[1/nSegments:1/nSegments:1-1/nSegments];
            points = 4*points-2;
            Melements=markElements(Mesh,[points]);
            New_Mesh=Mesh;
            New_Mesh = init_LEB(New_Mesh);
            New_Mesh = refine_LEB(New_Mesh,Melements);
            % Update mesh data structure
            New_Mesh = add_Edges(New_Mesh);
            New_Mesh = rmfield(New_Mesh,'BdFlags');
            New_Mesh= setBndFlags(New_Mesh);
            %plot_Mesh(New_Mesh)
    end
    
    New_Mesh=add_Edge2Elem(New_Mesh); 
    New_Mesh = add_DGData(New_Mesh);
    % topological derivatives
    if onlyStandard
        TopGrad=assemMat_TopGrad(New_Mesh);   % topological Gradient
    end
    TopRot=assemMat_TopRot(New_Mesh);      % topological Rotation

    % stiffness matrices
    M_W1F = assemMat_W1F(New_Mesh,@MASS_W1F,MU_HANDLE, P7O6());  % Mass Matrix Whitney-1-Forms
    C = assemMat_W1F(New_Mesh,@STIMA_Curl_W1F,MU_HANDLE, P7O6());  % Mass Matrix Whitney-1-Forms
    
    % boundary data on inflow boundary and righthandside
    c1=assemCochain_1f(New_Mesh,d.SOL1_Handle,gauleg(0,1,10));     % -v x curl u righthand side
    c2=assemCochain_1f(New_Mesh,d.SOL2_Handle,gauleg(0,1,10));     % grad(v.u) righthand side
    c3=assemCochain_1f(New_Mesh,d.SOL3_Handle,gauleg(0,1,10));     % curl curl rightand side
    c4=assemCochain_1f(New_Mesh,d.U_EX_Handle,gauleg(0,1,10));      % u righthand side
    
    Lie1 = assemMat_W1F(New_Mesh,@STIMA_ContrRot,d.V_Handle,P7O6());
    Lie2 =  assemMat_Inn_LieW1F(New_Mesh,@STIMA_Inn_LieW1F,d.V_Handle,QuadRule2D);
    %Lie2bnd1 =  assemMat_Bnd_LieW1F(New_Mesh,[-1, -2],@STIMA_Bnd_LieW1Fnormal,d.V_Handle,gauleg(0,1,10));
    %Lie2bnd2 =  assemMat_Bnd_LieW1F(New_Mesh,[-1, -2],@STIMA_Bnd_LieW1Ftangential,d.V_Handle,gauleg(0,1,5));
    %Lie2bnd=Lie2bnd1;
    Lie2bnd = assemMat_Bnd_LieW1F(New_Mesh,[-2],@STIMA_Bnd_LieW1F,d.V_Handle,gauleg(0,1,5));
    Lie2pen = assemMat_Inn_LieW1F(New_Mesh,@STIMA_Inn_LiePennW1F,d.V_Handle,QuadRule2D,0.5);
     % contraction
     if onlyStandard
         ContrOne=assemMat_Contr1f(New_Mesh,d.V_Handle);  % contraction of one forms
         V=d.V_Handle(New_Mesh.Coordinates);
         ContrTwo=assemMat_Contr2f(New_Mesh,V);   % contraction of two form
         A=M_W1F*ContrTwo*TopRot;               % -v x curl u geom.
         B=M_W1F*TopGrad*ContrOne;             % grad(v.u) geom.
     end
    
    % exact solution
    cU=assemCochain_1f(New_Mesh,d.U_EX_Handle,gauleg(0,1,10));
    cU_up= cU;
    cU_st = cU;
    cU_pst= cU;
        
    % calculate solutions
    % free DOFS
    nonInFlow=find(New_Mesh.BdFlags ~= -1);
    nonOutFlow=find(New_Mesh.BdFlags ~= -2);
    FreeDofs=nonInFlow;
    %FreeDofs=intersect( find(New_Mesh.BdFlags ~=-1), find(New_Mesh.BdFlags ~=-2));
    
    % standard scheme
    S_st=e_m*M_W1F+e_c*C+e_a*Lie1+e_b*Lie2-Lie2bnd;
    cf = M_W1F*(e_m*c4+e_a*c1+e_c*c3+e_b*c2);
    cfbnd = zeros(size(cf));
    cfbnd = assemLoad_Bnd_LieW1F(New_Mesh,[-1],cfbnd,gauleg(0,1,10),d.U_EX_Handle, d.V_Handle);
    cU_st(FreeDofs) = 0;
    Lf_st = cf+cfbnd-S_st*cU_st;
    % solve system
    cU_st(FreeDofs) = S_st(FreeDofs,FreeDofs)\Lf_st(FreeDofs);
    
    if onlyStandard
        % penalized standard scheme
        S_pst=e_m*M_W1F+e_c*C+e_a*Lie1+e_b*Lie2+e_b*Lie2pen-Lie2bnd;
        cf = M_W1F*(e_m*c4+e_a*c1+e_c*c3+e_b*c2);
        cfbnd = zeros(size(cf));
        cfbnd = assemLoad_Bnd_LieW1F(New_Mesh,[-1],cfbnd,gauleg(0,1,10),d.U_EX_Handle, d.V_Handle);
        %cfbnd = assemLoad_Bnd_LieW1Ftangential(New_Mesh,[-2],cfbnd,gauleg(0,1,10),d.U_EX_Handle, d.V_Handle);
        cU_pst(FreeDofs) = 0;
        Lf_pst = cf+cfbnd-S_pst*cU_pst;
        
        % solve system
        cU_pst(FreeDofs) = S_pst(FreeDofs,FreeDofs)\Lf_pst(FreeDofs);
    end
    
    % upwind scheme
    if onlyStandard
        S_up =e_m*M_W1F+e_c*C+e_a*A+e_b*B;%-e_bfB*Lie2bnd;
        %S_up =e_m*M_W1F+e_c*C+e_a*Lie1+e_b*B;%-e_bfB*Lie2bnd;
        cVUBd=assemBndCochain_0f(New_Mesh,[-2 -3],d.V_U_EX_Handle); 
        cf = M_W1F*(e_m*c4+e_a*c1+e_c*c3+e_b*c2-TopGrad*cVUBd);
        cU_up(FreeDofs) = 0;
        Lf_up = cf-S_up*cU_up;
        % solve system
        cU_up(FreeDofs)=(S_up(FreeDofs,FreeDofs))\Lf_up(FreeDofs);
    end
    %plot_Norm_W1F(cU_pst,New_Mesh); colorbar; title('U Patches');
    %  plot_P0(TopRot*cU_st,New_Mesh); colorbar; title('U Patches');
    if onlyStandard
        L2err1(i) = L2Err_W1F(New_Mesh,cU_pst,QuadRule,d.U_EX_Handle);
        Hcurlerr1(i) = HCurlSErr_W1F(New_Mesh,cU_pst,QuadRule,d.CURL_U_EX_Handle);
        IDerr1(i) = HveloCurlSErr_W1F(New_Mesh,cU_pst,QuadRule,d.SOL1_Handle,d.V_Handle);
        PIDerr1(i) = HPveloCurlSErr_W1F(New_Mesh,ContrTwo*TopRot*cU_pst,QuadRule,d.SOL1_Handle);
        L2err2(i) = L2Err_W1F(New_Mesh,cU_up,QuadRule,d.U_EX_Handle);
        Hcurlerr2(i) = HCurlSErr_W1F(New_Mesh,cU_up,QuadRule,d.CURL_U_EX_Handle);
        IDerr2(i) = HveloCurlSErr_W1F(New_Mesh,cU_up,QuadRule,d.SOL1_Handle,d.V_Handle);
        PIDerr2(i) = HPveloCurlSErr_W1F(New_Mesh,ContrTwo*TopRot*cU_up,QuadRule,d.SOL1_Handle);
    end
    L2err3(i) = L2Err_W1F(New_Mesh,cU_st,QuadRule,d.U_EX_Handle)
    Hcurlerr3(i) = HCurlSErr_W1F(New_Mesh,cU_st,QuadRule,d.CURL_U_EX_Handle);
    IDerr3(i) = HveloCurlSErr_W1F(New_Mesh,cU_st,QuadRule,d.SOL1_Handle,d.V_Handle);
    PIDerr3(i) = HPveloCurlSErr_W1F(New_Mesh,ContrTwo*TopRot*cU_st,QuadRule,d.SOL1_Handle);
    if onlyStandard
        JPen1(i) = sqrt(cU_pst'*Lie2pen*cU_pst); 
        JPen2(i) = sqrt(cU_up'*Lie2pen*cU_up); 
    end
    JPen3(i) = sqrt(cU_st'*Lie2pen*cU_st);
    
 %   [err1 err2 err3 err4 err5 err6 err7 err8 err9 JPen1 JPen2 JPen3]

    h(i) = get_MeshWidth(New_Mesh);
    Dofs(i) = size(New_Mesh.Edges,1);
    
end

%[diff(log(err1))./diff(log(h)) diff(log(err2))./diff(log(h)) diff(log(err3))./diff(log(h)) diff(log(err4))./diff(log(h)) diff(log(err5))./diff(log(h)) diff(log(err6))./diff(log(h))]
%clear all

fig = figure('Name','L^2-error');
plot(h,L2err3(1)/h(1)*h,'k--',h,L2err3(1)/h(1)^(1/2)*h.^(1/2),'k--',h,L2err3,'r-',h,L2err1,'b-',h,L2err2,'g-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf L^2-Error}');
legend('h','h^{1/2}','unstabelized','stabelized','Ext/Contr','Location','SouthEast')

fig = figure('Name','Penalty');
plot(h,JPen1(1)/h(1)*h,'k--',h,JPen1(1)/h(1)^(1/2)*h.^(1/2),'k--',h,JPen3,'r-',h,JPen1,'b-',h,JPen2,'g-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf facet-penalty-term}');
legend('h','h^{1/2}','unstabelized','stabelized','Ext/Contr','Location','SouthEast')

fig = figure('Name','Hcurl-error');
plot(h,Hcurlerr1(1)/h(1)*h,'k--',h,Hcurlerr1(1)/h(1)^(1/2)*h.^(1/2),'k--',h,Hcurlerr3,'r-',h,Hcurlerr1,'b-',h,Hcurlerr2,'g-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf semi Hcurl-Error}');
legend('h','h^{1/2}','unstabelized','stabelized','Ext/Contr','Location','SouthEast')

fig = figure('Name','veloHcurl-error');
plot(h,IDerr3(1)/h(1)*h,'k--',h,IDerr3(1)/h(1)^(1/2)*h.^(1/2),'k--',h,IDerr3,'r-',h,IDerr1,'b-',h,IDerr2,'g-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf veloHcurl-error}');
legend('h','h^{1/2}','unstabelized','stabelized','Ext/Contr','Location','SouthEast')

fig = figure('Name','veloHcurl-error');
plot(h,PIDerr3(1)/h(1)*h,'k--',h,PIDerr3(1)/h(1)^(1/2)*h.^(1/2),'k--',h,PIDerr3,'r-',h,PIDerr1,'b-',h,PIDerr2,'g-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf veloHcurl-error}');
legend('h','h^{1/2}','unstabelized','stabelized','Ext/Contr','Location','SouthEast')

toc
error = [L2err1 L2err2 L2err3]
slope = [diff(log(L2err1))./diff(log(h)) diff(log(L2err2))./diff(log(h)) diff(log(L2err3))./diff(log(h))]
save([filename,'.txt']','slope','h','error','-ASCII')

function MarkedElements=markElements(Mesh,ypoints)

m=length(ypoints);
nElements= size(Mesh.Elements,1);
MarkedElements=[];

for i=1:nElements
    x= 1/3*sum(Mesh.Coordinates(Mesh.Elements(i,:),1));
    y = 1/3*sum(Mesh.Coordinates(Mesh.Elements(i,:),2));
    if min(abs(ypoints-y+x))<eps;
        MarkedElements=[MarkedElements,i];
    end
end

function Mesh=setBndFlags(Mesh)

Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -2;
for i = Loc'
    if (1/2*sum(Mesh.Coordinates(Mesh.Edges(i,:),2))==-1 || 1/2*sum(Mesh.Coordinates(Mesh.Edges(i,:),1))==-1)
        Mesh.BdFlags(i) = -1;
    end
end