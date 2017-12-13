function main_SemiLagrangianSchemes(CFL,NREFS)
%  Semi-Lagrange- version for MHD

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%close all;

clear Mesh;
filename= ['./results/example1test_',num2str(CFL)];
%clear all;
NREFS_init = 0;     % Number of uniform red refinements
% NREFS = 3;
JIG =1;

MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;

%CFL =0.8; %
T0 = 0;
T1= 0.4;

mw = zeros(NREFS,1);
steps=zeros(NREFS,1);

c_CC =0;
c_Dt = 1;
c_ID = 1;
c_DI = c_ID;
    
T_Handle=@(t,varargin) 1;%cos(2*pi*t);              %1/(1+t);%1;%
DT_Handle=@(t,varargin) 0;%-2*pi*sin(2*pi*t);  %-1/(1+t)^2;%

% H1=@(x)sin(pi.*x(:,1)).*(1-x(:,2));
% D1H1=@(x)pi*(1-x(:,2)).*cos(pi.*x(:,1));
% D2H1=@(x)-sin(pi.*x(:,1));
% D11H1=@(x)-pi^2*(1-x(:,2)).*sin(pi.*x(:,1));
% D12H1=@(x)-pi.*cos(pi.*x(:,1));
% D22H1=@(x)0*x(:,1);
H1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1H1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2H1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
D11H1=@(x)-pi^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D12H1=@(x)pi^2*cos(pi.*x(:,2)).*cos(pi.*x(:,1));
D22H1=@(x)-pi^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));

H2=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
D1H2=@(x)(1-x(:,2).^2).*(-2*x(:,1));
D2H2=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
D11H2=@(x)(1-x(:,2).^2).*(-2);
D12H2=@(x)(-2*x(:,2)).*(-2*x(:,1));
D22H2=@(x)(-2).*(1-x(:,1).^2);

v1=0.66;
v2=1;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
%DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
%D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
%D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1DIR2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2DIR2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
%  DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
%  D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
%  D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
%  DIR2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
%  D1DIR2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
%  D2DIR2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
 
%  DIR1 = @(x)  (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
%  D1DIR1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
%  D2DIR1=@(x)  (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
%  DIR2=@(x) -1 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
%  D1DIR2=@(x) -1 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
%  D2DIR2=@(x) 4 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];
JvDirT_Handle=@(x,t)[D1DIR1(x) D2DIR1(x) D1DIR2(x) D2DIR2(x)];

H_Handle=@(x,flag,t,varargin)...
    T_Handle(t).*[H1(x) H2(x)];

% - curl h
CH_Handle=@(x,flag,t,varargin)T_Handle(t).*...
    D2H1(x) +D1H2(x);

% curl curl h
CCH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-D22H1(x) + D12H2(x) ...
    D12H1(x)-D11H2(x)];

% d_t h
DtH_HANDLE=@(x,flag,t,varargin)...
    DT_Handle(t).*[ H1(x) H2(x)];

%  -v x curl h
ID1H_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-DIR2(x).*(D1H2(x)-D2H1(x)) ...
       DIR1(x).*(D1H2(x)-D2H1(x))];

% grad v h
DI1H_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1DIR1(x).*H1(x)+DIR1(x).*D1H1(x)+D1DIR2(x).*H2(x)+DIR2(x).*D1H2(x) ...
    D2DIR1(x).*H1(x)+DIR1(x).*D2H1(x)+D2DIR2(x).*H2(x)+DIR2(x).*D2H2(x)];

% D h v
DHDIR_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1H1(x).*DIR1(x)+D2H1(x).*DIR2(x) ...
    D1H2(x).*DIR1(x)+D2H2(x).*DIR2(x)];

% D v h
DDIRH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1DIR1(x).*H1(x)+D2DIR1(x).*H2(x) ...
    D1DIR2(x).*H1(x)+D2DIR2(x).*H2(x)];

%  h x curl v
HDDIR_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [H2(x).*(D1DIR2(x)-D2DIR1(x)) ...
    -H1(x).*(D1DIR2(x)-D2DIR1(x))];

% v * div h
ID2H_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [DIR1(x).*(D1H1(x)+D2H2(x)) ...
    DIR2(x).*(D1H1(x)+D2H2(x))];

% -curl(v x h)
DI2H_HANDLE=@(x,flag,t,vvarargin)T_Handle(t).*...
    [-D2DIR1(x).*H2(x)-DIR1(x).*D2H2(x)+D2DIR2(x).*H1(x)+DIR2(x).*D2H1(x) ...
    +D1DIR1(x).*H2(x)+DIR1(x).*D1H2(x)-D1DIR2(x).*H1(x)-DIR2(x).*D1H1(x)];


%F_HANDLE=@(x,flag,t,varargin)...
%    c_CC * CCH_HANDLE(x,flag,t) + ...
%   c_Dt * DtH_HANDLE(x,flag,t) + ...
%    c_ID * ID1H_HANDLE(x,flag,t) + ...
%    c_DI * DI1H_HANDLE(x,flag,t);

F_HANDLE=@(x,flag,t,varargin)...
    c_CC * CCH_HANDLE(x,flag,t) + ...
    c_Dt * DtH_HANDLE(x,flag,t) + ...
    HDDIR_HANDLE(x,flag,t) + ...
    DDIRH_HANDLE(x,flag,t) + ...
    DHDIR_HANDLE(x,flag,t);

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

    h = CFL*mw(j)^(2/3)/norm([v1,v2]);
    nsteps = ceil((T1-T0)/h);
    h = h-(T0+nsteps*h-T1)/nsteps
    steps(j) = nsteps;

    NewMesh = init_LEB(NewMesh);
    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);
    NewMesh = add_DGData(NewMesh);
    
    % Assemble Curl-curl matrix, MASS matrix and load vector
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());  
    M2nd = assemMat_W1F2nd(Mesh,@MASS_W1F2nd,@(x,varargin)1, P7O6());
    
    % discrete pullback
    % Interpolation
    pbVm = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,-h);
    P_i = assemMat_SemiLag_W1F(NewMesh, pbVm);
    pbVp = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,h);

    % Quadrature
    pbBm = trace_bcenters(NewMesh, Dir_Handle,JvDir_Handle,-h);
    P_qbary = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbBm);
    
    P_qbary2nd = assemMat_SemiLagQuad_W1F_bary2nd(NewMesh, pbBm);

    % patches
    defMesh = NewMesh;
    defMesh.Coordinates = pbVp(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbVp);
    P_q = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    
    % patches 2nd
    P_pat2nd = assemMat_SemiLagQuad_W1F_patches2nd(NewMesh, defMesh,intersec);
    %P_pat2nd = assemMat_SemiLagQuad_W1F_patches2ndQuad(NewMesh, defMesh,intersec,P7O6());
        
    %Eulerian
    % patches
    Lie_p = assemMat_W1F(NewMesh,@STIMA_ContrRot,Dir_Handle,P7O6());
    Lie_p = Lie_p + assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LieW1F,Dir_Handle,gauleg(0,1,5));
    Lie_p = Lie_p + assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LiePennW1F,Dir_Handle,gauleg(0,1,5),0.5);
    % interpolation
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);      % topological Rotation
    ContrOne=assemMat_Contr1f(NewMesh,Dir_Handle);  % contraction of one forms
    V=Dir_Handle(NewMesh.Coordinates);
    ContrTwo=assemMat_Contr2f(NewMesh,V);   % contraction of two form
    Lie_i=M*ContrTwo*TopRot+M*TopGrad*ContrOne;

    % timestepping

    %[interpolation quadrature patches]
    H_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
%    H_test = assemLieAdCochain_1f(NewMesh,H_Handle,Dir_Handle,gauleg(0,1,1),0);
    H_exIP_old = H_old;
    H_exL2_old = H_old;
    H_old=[H_old H_old H_old H_old H_old];
    
    H2nd_old = M2nd\assemLoad_W1F2nd(Mesh,P7O6(),H_Handle,0);
    H2ndb_old = M2nd\assemLoad_W1F2nd(Mesh,P7O6(),H_Handle,0);
        
    time(1,j) = 0;
    L2ErrH(1,j,1) = L2Err_W1F(NewMesh,H_old(:,1),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,2) = L2Err_W1F(NewMesh,H_old(:,2),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,3) = L2Err_W1F(NewMesh,H_old(:,3),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,4) = L2Err_W1F(NewMesh,H_old(:,4),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,5) = L2Err_W1F(NewMesh,H_old(:,5),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,6) = L2Err_W1F2nd(NewMesh,H2nd_old,P7O6(),H_Handle,0,0);
    L2ErrH(1,j,7) = L2Err_W1F2nd(NewMesh,H2ndb_old,P7O6(),H_Handle,0,0);
    L2ErrH(1,j,8) = L2Err_W1F(NewMesh,H_exIP_old,P7O6(),H_Handle,0,0);
    L2ErrH(1,j,9) = L2Err_W1F(NewMesh,H_exL2_old,P7O6(),H_Handle,0,0);
    
    
    for i = 1:nsteps
        [i nsteps]
        time(i+1,j) = h+time(i,j);

        %Dirichlet data
        [H_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,5),(i)*h);
        H_new=[H_Dir H_Dir H_Dir H_Dir H_Dir];

        %system matrix
        A = c_Dt*M; 
        
        %righthand side
        Lf = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),i*h);
        L=[Lf Lf Lf Lf Lf];
        L = h*M*L+[M*P_i*H_old(:,1) P_qbary*H_old(:,2) P_q*H_old(:,3) (M-h*Lie_p)*H_old(:,4) (M-h*Lie_i)*H_old(:,5)]-A*H_new;
       
        %solve system
        H_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\L(FreeDofs,:);
       
        %Update vectors
        H_old=H_new;
       
        L2ErrH(i+1,j,1) = L2Err_W1F(NewMesh,H_old(:,1),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,2) = L2Err_W1F(NewMesh,H_old(:,2),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,3) = L2Err_W1F(NewMesh,H_old(:,3),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,4) = L2Err_W1F(NewMesh,H_old(:,4),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,5) = L2Err_W1F(NewMesh,H_old(:,5),P7O6(),H_Handle,0,i*h);
       
        % 2nd
        nEdges = size(NewMesh.Edges(),1);
        L2nd = assemLoad_W1F2nd(NewMesh,P7O6(),F_HANDLE,i*h);
        H2nd_new = zeros(2*nEdges,1);
        H2ndb_new = zeros(2*nEdges,1);
        H2nd_new = M2nd\(h*L2nd+P_pat2nd*H2nd_old);
        H2ndb_new = M2nd\(h*L2nd+P_qbary2nd*H2ndb_old);
        
        %Update vectors
        H2nd_old = H2nd_new;
        H2ndb_old = H2ndb_new;
        
        L2ErrH(i+1,j,6) = L2Err_W1F2nd(NewMesh,H2nd_old,P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,7) = L2Err_W1F2nd(NewMesh,H2ndb_old,P7O6(),H_Handle,0,i*h);
        
        % explizit
        [H_exIP_new,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,5),(i)*h);
        Lf = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),i*h);
        A= (M+h*Lie_i);
        L = h*M*Lf+M*H_exIP_old-A*H_exIP_new;
        H_exIP_new(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        L2ErrH(i+1,j,9) = L2Err_W1F(NewMesh,H_exIP_new,P7O6(),H_Handle,0,i*h);
        H_exIP_old = H_exIP_new;
        
        [H_exL2_new,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,5),(i)*h);
        Lf = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),i*h);
        A= (M+h*Lie_p);
        L = h* M*Lf+M*H_exL2_old-A*H_exL2_new;
        H_exL2_new(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        L2ErrH(i+1,j,8) = L2Err_W1F(NewMesh,H_exL2_new,P7O6(),H_Handle,0,i*h);
        H_exL2_old = H_exL2_new;
        
   end % timestep

end%  Mesh refinement

for j=1:NREFS
   j
   stopTimeL2ErrH(j,:)=L2ErrH(steps(j)+1,j,:)
end

figure;
hold on;
plot(mw,stopTimeL2ErrH(:,1),'rx-',...
    mw,stopTimeL2ErrH(:,2),'bx-',...
    mw,stopTimeL2ErrH(:,3),'gx-',...
    mw,stopTimeL2ErrH(:,6),'gx--',...
     mw,stopTimeL2ErrH(:,7),'bx--','Linewidth',2,'MarkerSize',8);
grid('on');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf L^2-error}','FontSize',12,'FontWeight','bold');
legend('Interpolation','Quadrature','Patch','Patch2nd','Quadrature2nd','Location','Northwest');
saveas(gcf,[filename,'error.fig'])
print('-depsc' ,[filename,'error.eps'])
hold off;

slope=[diff(log(stopTimeL2ErrH(:,1)))./diff(log(mw)) ...
    diff(log(stopTimeL2ErrH(:,2)))./diff(log(mw)) ...
    diff(log(stopTimeL2ErrH(:,3)))./diff(log(mw)) ...
    diff(log(stopTimeL2ErrH(:,4)))./diff(log(mw)) ...
    diff(log(stopTimeL2ErrH(:,5)))./diff(log(mw)) ...
    diff(log(stopTimeL2ErrH(:,6)))./diff(log(mw))...
    diff(log(stopTimeL2ErrH(:,7)))./diff(log(mw)) ...
    diff(log(stopTimeL2ErrH(:,8)))./diff(log(mw))...
    diff(log(stopTimeL2ErrH(:,9)))./diff(log(mw))]

figure;
hold on;
plot(mw,stopTimeL2ErrH(:,4),'rx-',...
    mw,stopTimeL2ErrH(:,5),'bx-',...
    mw,stopTimeL2ErrH(:,8),'rx--',...
    mw,stopTimeL2ErrH(:,9),'bx--','Linewidth',2,'MarkerSize',8);
grid('on');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf L^2-error}','FontSize',12,'FontWeight','bold');
legend('Patch explicit','Interpolation explicit','Patch implicit','Interpolation implicit','Location','Northwest');

save([filename,'.mat'],'slope','mw','stopTimeL2ErrH','CFL','H1','H2','DIR1','DIR2')
save([filename,'.txt']','slope','mw','stopTimeL2ErrH','CFL','-ASCII')
clear all

return