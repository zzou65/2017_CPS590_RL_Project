%function [time,mass,lm]=SemiLagrangeDirect(CFL)
function main_SemiLagrangianPureHeli(CFL,NREFS)
%  Semi-Lagrange- version for MHD

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%close all;
filename= ['./results/example3_',num2str(CFL)];
clear Mesh;
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

lambda=0.25;
H1=@(x)exp(-1/(2*lambda^2)*((x(:,1)+0.25).^2+(x(:,2)+0.25).^2));
H2=@(x)exp(-1/(2*lambda^2)*((x(:,1)+0.25).^2+(x(:,2)+0.25).^2));

H_Handle=@(x,flag,t,varargin)[H1(x) H2(x)];
F_HANDLE=@(x,flag,t,varargin)zeros(size(x));

v1=0.66;
v2=1;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];
JvDirT_Handle=@(x,t)[D1DIR1(x) D2DIR1(x) D1DIR2(x) D2DIR2(x)];

v1=-0.66;
v2=-1;
mDIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
mD1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
mD2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
mDIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
mD1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
mD2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

mDir_Handle=@(x,t)[mDIR1(x), mDIR2(x)];
mJvDir_Handle=@(x,t)[mD1DIR1(x) mD1DIR2(x) mD2DIR1(x) mD2DIR2(x)];
mJvDirT_Handle=@(x,t)[mD1DIR1(x) mD2DIR1(x) mD1DIR2(x) mD2DIR2(x)];

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
    nsteps = ceil((T1-T0)/(h));
    h = h-(T0+nsteps*h-T1)/nsteps
    steps(j) = 2*nsteps;

    NewMesh = init_LEB(NewMesh);
    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);
    NewMesh = add_DGData(NewMesh);
    
    % Assemble Curl-curl matrix, MASS matrix and load vector
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());  
   
    % positiv velocity
    
    % discrete pullback
    % Interpolation
    pbVm = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,-h);
    P_i = assemMat_SemiLag_W1F(NewMesh, pbVm);
    pbVp = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,h);

    % Quadrature
    pbBm = trace_bcenters(NewMesh, Dir_Handle,JvDir_Handle,-h);
    P_qbary = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbBm);
    
    % patches
    defMesh = NewMesh;
    defMesh.Coordinates = pbVp(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbVp);
    P_q = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
            
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
    H_ref=H_old;
    H_old=[H_old H_old H_old H_old H_old];
    %plot_Norm_W1F(H_old(:,1),NewMesh);
         
    time(1,j) = 0;
    L2ErrH(1,j,1) = L2Err_W1F(NewMesh,H_old(:,1),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,2) = L2Err_W1F(NewMesh,H_old(:,2),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,3) = L2Err_W1F(NewMesh,H_old(:,3),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,4) = L2Err_W1F(NewMesh,H_old(:,4),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,5) = L2Err_W1F(NewMesh,H_old(:,5),P7O6(),H_Handle,0,0);
    
    mass(j,1) = L2Err_W1F(NewMesh,H_old(:,1),P7O6(),F_HANDLE,0,0);
    mass(j,2) = L2Err_W1F(NewMesh,H_old(:,2),P7O6(),F_HANDLE,0,0);
    mass(j,3) = L2Err_W1F(NewMesh,H_old(:,3),P7O6(),F_HANDLE,0,0);
    mass(j,4) = L2Err_W1F(NewMesh,H_old(:,4),P7O6(),F_HANDLE,0,0);
    mass(j,5) = L2Err_W1F(NewMesh,H_old(:,5),P7O6(),F_HANDLE,0,0);
    
    heli(1,j,1) = Helicity2D_W1F(NewMesh,H_old(:,1),P7O6(),Dir_Handle,0);
    heli(1,j,2) = Helicity2D_W1F(NewMesh,H_old(:,2),P7O6(),Dir_Handle,0);
    heli(1,j,3) = Helicity2D_W1F(NewMesh,H_old(:,3),P7O6(),Dir_Handle,0);
    heli(1,j,4) = Helicity2D_W1F(NewMesh,H_old(:,4),P7O6(),Dir_Handle,0);
    heli(1,j,5) = Helicity2D_W1F(NewMesh,H_old(:,5),P7O6(),Dir_Handle,0); 
    
    for i = 1:nsteps
        [i 2*nsteps]
        time(i+1,j) = h+time(i,j);

        %Dirichlet data
        [H_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,5),(i)*h);
        H_new=[H_Dir H_Dir H_Dir H_Dir H_Dir];

        %system matrix
        A = M; 
        
        %righthand side
        Lf = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),i*h);
        L=[Lf Lf Lf Lf Lf];
        L = h*M*L+[M*P_i*H_old(:,1) P_qbary*H_old(:,2) P_q*H_old(:,3) (M-Lie_p)*H_old(:,4) (M-Lie_i)*H_old(:,5)]-A*H_new;
       
        %solve system
        H_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\L(FreeDofs,:);
       
        %Update vectors
        H_old=H_new;
        
        L2ErrH(i+1,j,1) = L2Err_W1F(NewMesh,H_old(:,1),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,2) = L2Err_W1F(NewMesh,H_old(:,2),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,3) = L2Err_W1F(NewMesh,H_old(:,3),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,4) = L2Err_W1F(NewMesh,H_old(:,4),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,5) = L2Err_W1F(NewMesh,H_old(:,5),P7O6(),H_Handle,0,i*h);
    
        heli(i+1,j,1) = Helicity2D_W1F(NewMesh,H_old(:,1),P7O6(),Dir_Handle,0);
        heli(i+1,j,2) = Helicity2D_W1F(NewMesh,H_old(:,2),P7O6(),Dir_Handle,0);
        heli(i+1,j,3) = Helicity2D_W1F(NewMesh,H_old(:,3),P7O6(),Dir_Handle,0);
        heli(i+1,j,4) = Helicity2D_W1F(NewMesh,H_old(:,4),P7O6(),Dir_Handle,0);
        heli(i+1,j,5) = Helicity2D_W1F(NewMesh,H_old(:,5),P7O6(),Dir_Handle,0);

   end % timestep
   
   % negative velocity
   
    % discrete pullback
    % Interpolation
    pbVm = trace_vertices_W(NewMesh,mDir_Handle,mJvDir_Handle,-h);
    P_i = assemMat_SemiLag_W1F(NewMesh, pbVm);
    pbVp = trace_vertices_W(NewMesh,mDir_Handle,mJvDir_Handle,h);

    % Quadrature
    pbBm = trace_bcenters(NewMesh, mDir_Handle,mJvDir_Handle,-h);
    P_qbary = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbBm);
    
    % patches
    defMesh = NewMesh;
    defMesh.Coordinates = pbVp(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbVp);
    P_q = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    
    %Eulerian
    % patches
    Lie_p = assemMat_W1F(NewMesh,@STIMA_ContrRot,mDir_Handle,P7O6());
    Lie_p = Lie_p + assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LieW1F,mDir_Handle,gauleg(0,1,5));
    Lie_p = Lie_p + assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LiePennW1F,mDir_Handle,gauleg(0,1,5),0.5);
    % interpolation
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);      % topological Rotation
    ContrOne=assemMat_Contr1f(NewMesh,mDir_Handle);  % contraction of one forms
    V=mDir_Handle(NewMesh.Coordinates);
    ContrTwo=assemMat_Contr2f(NewMesh,V);   % contraction of two form
    Lie_i=M*ContrTwo*TopRot+M*TopGrad*ContrOne;

    % timestepping

    %[interpolation quadrature patches]
    time(1,j) = 0;
    L2ErrH(1+nsteps,j,1) = L2Err_W1F(NewMesh,H_old(:,1),P7O6(),H_Handle,0,0);
    L2ErrH(1+nsteps,j,2) = L2Err_W1F(NewMesh,H_old(:,2),P7O6(),H_Handle,0,0);
    L2ErrH(1+nsteps,j,3) = L2Err_W1F(NewMesh,H_old(:,3),P7O6(),H_Handle,0,0);
    L2ErrH(1+nsteps,j,4) = L2Err_W1F(NewMesh,H_old(:,4),P7O6(),H_Handle,0,0);
    L2ErrH(1+nsteps,j,5) = L2Err_W1F(NewMesh,H_old(:,5),P7O6(),H_Handle,0,0);

    heli(1+nsteps,j,1) = Helicity2D_W1F(NewMesh,H_old(:,1),P7O6(),Dir_Handle,0);
    heli(1+nsteps,j,2) = Helicity2D_W1F(NewMesh,H_old(:,2),P7O6(),Dir_Handle,0);
    heli(1+nsteps,j,3) = Helicity2D_W1F(NewMesh,H_old(:,3),P7O6(),Dir_Handle,0);
    heli(1+nsteps,j,4) = Helicity2D_W1F(NewMesh,H_old(:,4),P7O6(),Dir_Handle,0);
    heli(1+nsteps,j,5) = Helicity2D_W1F(NewMesh,H_old(:,5),P7O6(),Dir_Handle,0);

    for i = (1+nsteps):2*nsteps
        [i 2*nsteps]
        time(i+1,j) = h+time(i,j);

        %Dirichlet data
        [H_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,5),(i)*h);
        H_new=[H_Dir H_Dir H_Dir H_Dir H_Dir];

        %system matrix
        A = M; 
        
        %righthand side
        Lf = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),i*h);
        L=[Lf Lf Lf Lf Lf];
        L = h*M*L+[M*P_i*H_old(:,1) P_qbary*H_old(:,2) P_q*H_old(:,3) (M-Lie_p)*H_old(:,4) (M-Lie_i)*H_old(:,5)]-A*H_new;
       
        %solve system
        H_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\L(FreeDofs,:);
       
        %Update vectors
        H_old=H_new;
%         plot_Norm_W1F(H_old(:,1),NewMesh);
       
        L2ErrH(i+1,j,1) = L2Err_W1F(NewMesh,H_old(:,1),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,2) = L2Err_W1F(NewMesh,H_old(:,2),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,3) = L2Err_W1F(NewMesh,H_old(:,3),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,4) = L2Err_W1F(NewMesh,H_old(:,4),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,5) = L2Err_W1F(NewMesh,H_old(:,5),P7O6(),H_Handle,0,i*h);

        heli(i+1,j,1) = Helicity2D_W1F(NewMesh,H_old(:,1),P7O6(),Dir_Handle,0);
        heli(i+1,j,2) = Helicity2D_W1F(NewMesh,H_old(:,2),P7O6(),Dir_Handle,0);
        heli(i+1,j,3) = Helicity2D_W1F(NewMesh,H_old(:,3),P7O6(),Dir_Handle,0);
        heli(i+1,j,4) = Helicity2D_W1F(NewMesh,H_old(:,4),P7O6(),Dir_Handle,0);
        heli(i+1,j,5) = Helicity2D_W1F(NewMesh,H_old(:,5),P7O6(),Dir_Handle,0);

   end % timestep
   
    mass(j,1) = abs(mass(j,1)-L2Err_W1F(NewMesh,H_old(:,1),P7O6(),F_HANDLE,0,0))/mass(j,1);
    mass(j,2) = abs(mass(j,2)-L2Err_W1F(NewMesh,H_old(:,2),P7O6(),F_HANDLE,0,0))/mass(j,2);
    mass(j,3) = abs(mass(j,3)-L2Err_W1F(NewMesh,H_old(:,3),P7O6(),F_HANDLE,0,0))/mass(j,3);
    mass(j,4) = abs(mass(j,4)-L2Err_W1F(NewMesh,H_old(:,4),P7O6(),F_HANDLE,0,0))/mass(j,4);
    mass(j,5) = abs(mass(j,5)-L2Err_W1F(NewMesh,H_old(:,5),P7O6(),F_HANDLE,0,0))/mass(j,5);

    %plot_Norm_W1F(H_old(:,1)-H_ref,NewMesh);
end%  Mesh refinement

for j=1:NREFS
   stopTimeL2ErrH(j,:)=L2ErrH(steps(j)+1,j,:);
   stopTimeheli(j,:)=heli(steps(j)+1,j,:);
end

figure;
hold on;
plot(mw,stopTimeL2ErrH(:,1),'rx-',...
    mw,stopTimeL2ErrH(:,2),'bx-',...
    mw,stopTimeL2ErrH(:,3),'gx-',...
    mw,stopTimeL2ErrH(:,4),'gx--',...
     mw,stopTimeL2ErrH(:,5),'bx--','Linewidth',2,'MarkerSize',8);
grid('on');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf L^2-error}','FontSize',12,'FontWeight','bold');
legend('Interpolation','Quadrature','Patch','Eulp','Eulip','Location','Northwest');
saveas(gcf,[filename,'error.fig'])
print('-depsc' ,[filename,'error.eps'])

figure;
hold on;
plot(mw,stopTimeheli(:,1),'rx-',...
    mw,stopTimeheli(:,2),'bx-',...
    mw,stopTimeheli(:,3),'gx-',...
    mw,stopTimeheli(:,4),'gx--',...
     mw,stopTimeheli(:,5),'bx--','Linewidth',2,'MarkerSize',8);
grid('on');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf helicity}','FontSize',12,'FontWeight','bold');
legend('Interpolation','Quadrature','Patch','Eulp','Eulip','Location','Northwest');
saveas(gcf,[filename,'error.fig'])
print('-depsc' ,[filename,'error.eps'])

figure;
hold on;
plot(mw,mass(:,1),'rx-',...
    mw,mass(:,2),'bx-',...
    mw,mass(:,3),'gx-',...
    mw,mass(:,4),'gx--',...
     mw,mass(:,5),'bx--','Linewidth',2,'MarkerSize',8);
grid('on');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf diffusivity}','FontSize',12,'FontWeight','bold');
legend('Interpolation','Quadrature','Patch','Eulp','Eulip','Location','Northwest');

saveas(gcf,[filename,'mass.fig'])
print('-depsc' ,[filename,'mass.eps'])

slope=[diff(log(stopTimeL2ErrH(:,1)))./diff(log(mw)) ...
    diff(log(stopTimeL2ErrH(:,2)))./diff(log(mw)) ...
    diff(log(stopTimeL2ErrH(:,3)))./diff(log(mw)) ...
    diff(log(stopTimeL2ErrH(:,4)))./diff(log(mw)) ...
    diff(log(stopTimeL2ErrH(:,5)))./diff(log(mw))]

slope=[diff(log(stopTimeheli(:,1)))./diff(log(mw)) ...
    diff(log(stopTimeheli(:,2)))./diff(log(mw)) ...
    diff(log(stopTimeheli(:,3)))./diff(log(mw)) ...
    diff(log(stopTimeheli(:,4)))./diff(log(mw)) ...
    diff(log(stopTimeheli(:,5)))./diff(log(mw))]

save([filename,'.mat'],'slope','mw','stopTimeL2ErrH','mass','CFL','H1','H2','DIR1','DIR2')
save([filename,'.txt']','slope','mw','stopTimeL2ErrH','mass','CFL','-ASCII')
clear all;