%function [time,mass,lm]=SemiLagrangeDirect(CFL)
%  Semi-Lagrange- version for MHD

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%close all;
clear Mesh;
%clear all;
NREFS_init = 0;% Number of uniform red refinements
NREFS = 4;
JIG =2;
PLOT=0;
THETA = 0 ; %0 explicit
if PLOT
    ! rm ./vtk/h0*
end

%
MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;
%

CFL =0.4;
T0 = 0;
T1= 0.5;

mw = zeros(NREFS,1);
Err = [];
time = [];
div=[];
mass=[];
steps=zeros(NREFS,1);

%
c_CC = 1;
c_Dt = 1;
c_ID = 1;
c_DI = c_ID;

T_Handle=@(t,varargin) cos(2*pi*t);%1/(1+t);%1;%
DT_Handle=@(t,varargin) -2*pi*sin(2*pi*t);%-1/(1+t)^2;%

% H1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D1H1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
% D2H1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
% D11H1=@(x)-pi.^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D12H1=@(x)pi^2*cos(pi.*x(:,2)).*cos(pi.*x(:,1));
% D22H1=@(x)-pi^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% H2=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
% D1H2=@(x)(1-x(:,2).^2).*(-2*x(:,1));
% D2H2=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
% D11H2=@(x)(1-x(:,2).^2).*(-2);
% D12H2=@(x)(-2*x(:,2)).*(-2*x(:,1));
% D22H2=@(x)(-2).*(1-x(:,1).^2);

% H1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D1H1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
% D2H1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
% D11H1=@(x)-pi.^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D12H1=@(x)pi^2*cos(pi.*x(:,2)).*cos(pi.*x(:,1));
% D22H1=@(x)-pi^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% H2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D1H2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
% D2H2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
% D11H2=@(x)-pi.^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D12H2=@(x)pi^2*cos(pi.*x(:,2)).*cos(pi.*x(:,1));
% D22H2=@(x)-pi^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));

% H1 = @(x) 4 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
% D1H1=@(x) -16 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
% D2H1=@(x) 4 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
% D11H1=@(x) -16 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
% D12H1=@(x) -16 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
% D22H1=@(x) -24 * (1-x(:,1).^2).^2 .*x(:,2);
% H2=@(x) -4 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
% D1H2=@(x) -4 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
% D2H2=@(x) 16 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
% D11H2=@(x) 24 * x(:,1).*(1-x(:,2).^2).^2;
% D12H2=@(x) 16 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
% D22H2=@(x) 16 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);

% H1=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
% D1H1=@(x)(1-x(:,2).^2).*(-2*x(:,1));
% D11H1=@(x)(1-x(:,2).^2)*(-2);
% D12H1=@(x)(-2*x(:,2)).*(-2*x(:,1));
% D2H1=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
% D22H1=@(x)(-2).*(1-x(:,1).^2);
% H2=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
% D1H2=@(x)(1-x(:,2).^2).*(-2*x(:,1));
% D11H2=@(x)(1-x(:,2).^2)*(-2);
% D12H2=@(x)(-2*x(:,2)).*(-2*x(:,1));
% D2H2=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
% D22H2=@(x)(-2).*(1-x(:,1).^2);

v1=0.66;
v2=1;
% DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
% D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
% D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
% DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
% D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
% D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);
% 
% DIR1=@(x)sin(pi*x(:,1)).*(1-x(:,2).^2);
% D1DIR1=@(x)pi*cos(pi*x(:,1)).*(1-x(:,2).^2);
% D2DIR1=@(x)sin(pi*x(:,1)).*(-2*x(:,2));
% DIR2=@(x)sin(pi*x(:,2)).*(1-x(:,1).^2);
% D1DIR2=@(x)sin(pi*x(:,2)).*(-2*x(:,1));
% D2DIR2=@(x)pi*cos(pi*x(:,2)).*(1-x(:,1).^2);
 
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1DIR2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2DIR2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));

% DIR1=@(x)v1*(1-x(:,2).^2).*(x(:,1).^2-x(:,1).^4);
% D1DIR1=@(x)v1*(1-x(:,2).^2).*(2*x(:,1)-4*x(:,1).^3);
% D2DIR1=@(x)v1*(-2*x(:,2)).*(x(:,1).^2-x(:,1).^4);
% DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
% D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
% D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];

H_Handle=@(x,flag,t,varargin)...
    T_Handle(t).*[ H1(x) H2(x)];

% curl curl h
CCH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-D22H1(x) + D12H2(x) ...
    D12H1(x)-D11H2(x)];

% d_t h
DtH_HANDLE=@(x,flag,t,varargin)...
    DT_Handle(t).*[ H1(x) H2(x)];

%  -v x curl h
IDH_HANDLE=@(x,flag,t,varargin)1.*T_Handle(t).*...
    [-DIR2(x).*(D1H2(x)-D2H1(x)) ...
      DIR1(x).*(D1H2(x)-D2H1(x))];

% grad v h
DIH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1DIR1(x).*H1(x)+DIR1(x).*D1H1(x)+D1DIR2(x).*H2(x)+DIR2(x).*D1H2(x) ...
     D2DIR1(x).*H1(x)+DIR1(x).*D2H1(x)+D2DIR2(x).*H2(x)+DIR2(x).*D2H2(x)];

F_HANDLE=@(x,flag,t,varargin)...
    c_CC * CCH_HANDLE(x,flag,t) + ...
    c_Dt * DtH_HANDLE(x,flag,t) + ...
    c_ID * IDH_HANDLE(x,flag,t) + ...
    c_DI * DIH_HANDLE(x,flag,t);

% Dir_Handle=@(x,t)[(1-x(:,2).^2).*(1-x(:,1).^2).*(-2*x(:,2)).*(1-x(:,1).^2)...
%     (1-x(:,2).^2).*(1-x(:,1).^2).*(1-x(:,2).^2).*(2*x(:,1))];
% Dir_Handle=@(x,t)straight_flow(x);%circ_flow(x);%
% F_HANDLE=@(x,flag,t,varargin)zeros(size(x,1),2);
% H_Handle=@(x,flag,t,varargin)...
%     [pulse_2D(x) pulse_2D(x)];
         
% Load mesh from file
Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

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
    
    mw(j) = get_MeshWidth(NewMesh);

    h = CFL*mw(j)/norm([v1,v2]);
    nsteps = ceil((T1-T0)/h);
    h = h-(T0+nsteps*h-T1)/nsteps;
    steps(j) = nsteps;

    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);

    % plot_Mesh(NewMesh,'petas');

    nEdges= size(NewMesh.Edges,1);

    % Assemble Curl-curl matrix, MASS matrix and load vector
    C = assemMat_W1F(NewMesh,@STIMA_Curl_W1F,SIGMA_HANDLE,P7O6());
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());
    
    % topological derivatives
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);     % topological Rotation

    % timestepping
    H_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    
    %Edge Length
    EdgeLength=Edge_Length(NewMesh);
    
    time(1,j) = 0;
    Err(1,j) = L2Err_W1F(NewMesh,H_old,P7O6(),H_Handle,0,0);
    div(1,j) = norm(TopGrad'*M*H_old);
    tens_unit = norm(EdgeLength.*H_old);
    mass_unit = norm(H_old'*M*H_old);
    mass(1,j) = 1;
    tension(1,j) =1;
    
    % pullback of edges
    directions = Dir_Handle(NewMesh.Coordinates,i*h);
    pbV = trace_vertices(NewMesh,-h*directions);
    P = assemMat_SemiLag_W1F(NewMesh, pbV);
    %pbV = trace_bcenters(NewMesh, Dir_Handle,JvDir_Handle,-h);
    %P = assemMat_SemiLagQuad_W1F(NewMesh, pbV);                                         % vertex
    %pbB = trace_bcenters(NewMesh, Dir_Handle,JvDir_Handle,-h);
    %P = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbB);                                  % barycenter
    P = assemMat_SemiLagQuad_W1F_strang1(NewMesh, @(x)Dir_Handle(x,0),JvDir_Handle,-h);      % strang1
    
    for i = 1:nsteps
        [i nsteps]
        time(i+1,j) = h+time(i,j);
       
        [H_new,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,1),(i)*h);

        L = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),(i)*h);
        
        A = c_Dt*M + c_CC*h*C;
        %L = h*M*L+M*P*H_old-A*H_new;
        L = h*M*L+P*H_old-A*H_new;
        
        H_new(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        
        % Update vectors
        H_old=H_new;
        
        if PLOT
            filename=sprintf('./vtk/h%05d',i);
            plot_W1F_vtk(H_old,NewMesh,filename);
        end

        Err(i+1,j) = L2Err_W1F(NewMesh,H_old,P7O6(),H_Handle,0,i*h);
        [Dummy,FreeDofs] = assemDir_LFE(NewMesh,-1,@(x,varargin)zeros(size(x,1),1));
        div(i+1,j) = norm(TopGrad(:,FreeDofs)'*M*H_old);
        mass(i+1,j) = norm(H_old'*M*H_old)/mass_unit;
        tension(i+1,j) = norm(EdgeLength.*H_old)/tens_unit;

         %plot_Norm_W1F(H_old, NewMesh); colorbar;
    end % timestep
    %figure;
    %plot_W1F(H_old, NewMesh);
   

   
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

fig = figure('Name','weak divergence');
for j=1:NREFS
    plot(time(1:steps(j)+1,j),div(1:steps(j)+1,j));
    hold on;
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
ylabel('{\bf weak divergence}');
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

fig = figure('Name','tension conservation');
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