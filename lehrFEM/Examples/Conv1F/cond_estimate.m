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
JIG=2;                   % 
NREFS =6;             % number of refinements
QuadRule = P7O6();

% model parameters
isigma = 10.^-[7];             % diffusion constant
mu = 1                  %
supg1=1;                 %impact of supg modification
supg2=1;                 %impact of supg modification

% select test code
d=getData_EsHv();

% Initialize mesh
Mesh.Coordinates = d.Coordinates;
Mesh.Elements = d.Elements;
Mesh = add_Edges(Mesh);

Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = d.boundtype;
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);

h=zeros(NREFS,1);
condUp=zeros(NREFS,length(isigma));
condUp=zeros(NREFS,length(isigma));

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

    % Stiffness matrices,
    % upwind
    DD = TopRot'*MassTwo*TopRot;             % curl curl u geom.
    ID_up = M_W1F*ContrTwo*TopRot;          % -v x curl u geom.

    DI_up = M_W1F*TopGrad*ContrOne;      % grad(v.u) geom.
    %ID_st = GradContr;                               % grad(v.u) geom.
    %ID_up = GradContr;                              % grad(v.u) geom.

    % Direchlet boundary
    [Hv_up,FreeDofs] = assemDir_W1F(NewMesh,-1,d.Hv_Handle,gauleg(0,1,10),1);

    for j =1:length(isigma)
        %  stiffnes matrices
        A_up =(isigma(j)*DD-ID_up'+M_W1F);
        A_upc =(isigma(j)*DD-ID_up'-DI_up'+M_W1F);
        % A_supg =(A_st+D_supg);

        condUp(i,j) = condest(A_up(FreeDofs,FreeDofs));
        condUpC(i,j) = condest(A_upc(FreeDofs,FreeDofs));
        %condUp(i,j) = cond(full(A_up(FreeDofs,FreeDofs)));
        %condUpC(i,j) = cond(full(A_upc(FreeDofs,FreeDofs)));
    end
end

figure()
plot(h,h.^(-2)*condUp(1,end)/h(1),'o-',h,condUp(:,end),'o-',h,condUpC(:,end),'o-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('h');
ylabel('cond');
legend('h^{-2}','Up','UpC','Location','northwest');

figure()
plot(isigma,isigma.^(-1)*condUp(end,1)/isigma(1),'o-',isigma,condUp(end,:),'o-',isigma,condUpC(end,:),'o-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('\sigma^{-1}');
ylabel('cond');
legend('h^{-1}','Up','UpC','Location','northwest');

