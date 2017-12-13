%  Semi-Lagrange- version for vorticity
%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

! rm ./vtk/vort0*
! rm ./vtk/phi0*
! rm ./vtk/velo0*

%close all;
clear all;
clear Mesh;
NREFS_init = 2;% Number of uniform red refinements-
NREFS = 1;
JIG =1;
THETA=1;    % 0 explizit

%
MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;
%

CFL = 1;
T0 = 0;
T1= 4;

mw = zeros(NREFS,1);
Err = [];
time = [];
div=[];
steps=zeros(NREFS,1);

v1=1;
v2=1;

DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
Dir_Handle=@(x,t)[(1-x(:,2).^2).*(1-x(:,1).^2).*(-2*x(:,2)).*(1-x(:,1).^2)...
    (1-x(:,2).^2).*(1-x(:,1).^2).*(1-x(:,2).^2).*(2*x(:,1))];
F_HANDLE=@(x,flag,t,varargin)zeros(size(x,1),1);
H_Handle=@(x,flag,t,varargin) pulse_2D(x);

% Load mesh from file

Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure

Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);

% refine
for i=1:NREFS_init
    Mesh=refine_REG(Mesh);
end

for j=1:NREFS
     Mesh=refine_REG(Mesh);
     % Mesh preprocessing
     
    switch (JIG)
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
    
    % CFL number
    mw(j) = get_MeshWidth(NewMesh);

    h = CFL*mw(j)/norm([v1,v2]);
    nsteps = ceil((T1-T0)/h);
    h = h-(T0+nsteps*h-T1)/nsteps;
    steps(j) = nsteps;

    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);

    % plot_Mesh(NewMesh,'petas');

    nEdges= size(NewMesh.Edges,1);
    
    % Assemble  MASS matrix 
    M = assemMat_P0(NewMesh, @MASS_P0);
    
    A = assemMat_LFE(NewMesh, @STIMA_Lapl_LFE);
    
    M_LFE_P0 = assemMat_P1P0(NewMesh, @MASS_P1P0);
    
    % topological derivatives
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);     % topological Rotation
    
    % contraction operators
    %ContrOne=assemMat_Contr1f(NewMesh,Dir_Handle);  % contraction of one forms
    V=Dir_Handle(NewMesh.Coordinates);
    ContrTwo=assemMat_Contr2f(NewMesh,V);   % contraction of two forms

    %ID = ContrOne*TopGrad;              % -v x curl u geom.
    DI = TopRot*ContrTwo;             % grad(v.u) geom.

    S_left = M+h*THETA*M*DI;
    S_right = M-h*(1-THETA)*M*DI;
    
    % timestepping

    H_init = assemLoad_P0(NewMesh,P7O6(),H_Handle,0);
    H_old = M\H_init;
    L_old = assemLoad_P0(NewMesh, P7O6(), F_HANDLE,0);

    time(1,j) = 0;
    Err(1,j) = L2Err_PC(NewMesh,H_old,P7O6(),H_Handle,0,0);
    mass_unit = norm(H_old'*M*H_old);
    mass(1,j) = 1;
    
    filename=sprintf('./vtk/vort%05d',0);
    plot_P0_vtk(H_old,NewMesh,filename);

    [phi, FDP1]=assemDir_LFE(NewMesh,-1,@(x,varargin)zeros(size(x,1),1));
    phi(FDP1)=A(FDP1,FDP1)\(M_LFE_P0(FDP1,FDP1)*H_old(FDP1));
    filename=sprintf('./vtk/phi%05d',0);
    plot_LFE_vtk(phi,NewMesh,filename);
    
    filename=sprintf('./vtk/velo%05d',0);
    plot_W1F_vtk(TopGrad*phi,NewMesh,filename);
    
    for i = 1:nsteps
        i
        time(i+1,j) = h+time(i,j);
        
        L_new = assemLoad_P0(NewMesh, P7O6(), F_HANDLE,i*h);
       
        L = h*(1-THETA)*L_old + h*THETA*L_new+S_right*H_old;

        H_new = S_left\L;
        
        % Update vectors
        H_old = H_new;
        L_old = L_new;

        filename=sprintf('./vtk/vort%05d',i);
        plot_P0_vtk(H_old,NewMesh,filename);
    
        [phi, FDP1]=assemDir_LFE(NewMesh,-1,@(x,varargin)zeros(size(x,1),1));
        phi(FDP1)=A(FDP1,FDP1)\(M_LFE_P0(FDP1,FDP1)*H_old(FDP1));
        filename=sprintf('./vtk/phi%05d',i);
        plot_LFE_vtk(phi,NewMesh,filename);
        
        filename=sprintf('./vtk/velo%05d',i);
        plot_W1F_vtk(TopGrad*phi,NewMesh,filename);
        
        Err(i+1,j) = L2Err_PC(NewMesh,H_old,P7O6(),H_Handle,0,i*h);
        mass(i+1,j) = norm(H_old'*M*H_old)/mass_unit;

    end % timestep
    %     figure;
    %     plot_W1F(H_old, NewMesh);

   
end%  Mesh refinement

finErr=zeros(NREFS,1);
fig = figure('Name','Error material derivative');
leg=[];
for j=1:NREFS
    plot(time(1:steps(j)+1,j),Err(1:steps(j)+1,j));
    hold on;
    finErr(j) = Err(steps(j)+1,j);
    lm(j,:)=['h=',sprintf('%1.3f',mw(j)),...
        ' \Delta t=',sprintf('%1.3f',(T1-T0)/steps(j)),...
        ' CFL=',sprintf('%1.3f',norm([v1 v2])*(T1-T0)/(mw(j)*steps(j)))];
end
grid('on');
set(gca,'YScale','log');
Markers = '.ox+*sdv^<>ph';
noMarkers = length(Markers);
Colors = 'gbrcmy';
noColors = length(Colors);
H1c = findobj(gca,'Type','line');
linehandles = [H1c];
for K = 1:length(linehandles)
    set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
    set(linehandles(K),'Color',Colors(1+mod(K-1,noColors)));
end
legend(lm);
xlabel('{\bf time}');
ylabel('{\bf error}');
hold off;

fig = figure('Name','mass conservation');
for j=1:NREFS
    plot(time(1:steps(j)+1,j),mass(1:steps(j)+1,j));
    hold on;
end
grid('on');
%set(gca,'YScale','log');
Markers = '.ox+*sdv^<>ph';
noMarkers = length(Markers);
Colors = 'gbrcmy';
noColors = length(Colors);
H1c = findobj(gca,'Type','line');
linehandles = [H1c];
for K = 1:length(linehandles)
    set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
    set(linehandles(K),'Color',Colors(1+mod(K-1,noColors)));
end
legend(lm);
xlabel('{\bf time}');
ylabel('{\bf mass}');
hold off;

figure;
plot(mw,finErr);
grid('on');
set(gca,'YScale','log','XScale','log');
p = polyfit(log(mw),log(finErr),1);
add_Slope(gca,'East',p(1));