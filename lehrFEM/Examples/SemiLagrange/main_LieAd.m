function main_LieAd
%  Pure advection, via time-integration, interpolation, stabelized standard scheme.
%  !!! We do not see convergence, since the analytical is not smooth across
%  the trajectroy of (-1,-1)

%   Copyright 2010 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%close all;

clear Mesh;
%clear all;
NREFS_init = 0;     % Number of uniform red refinements
NREFS =6;
JIG =1;

MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;

mw = zeros(NREFS,1);
L2ErrH=zeros(NREFS,3);
JPen=zeros(NREFS,3);

%H1=@(x) sin(pi/2.*x(:,2)).*sin(pi/2.*x(:,1));
%H2=@(x) cos(pi/2.*x(:,2)).*cos(pi/2.*x(:,1));

H1=@(x) sigmoid(x(:,1),-0.25).*(sin(pi/2.*x(:,2)).*sin(pi/2.*x(:,1)))+...
   sigmoid(x(:,2),-0.25).*(sin(pi/2.*x(:,2)).*sin(pi/2.*x(:,1)));
H2=@(x) sigmoid(x(:,1),-0.25).*(cos(pi/2.*x(:,2)).*cos(pi/2.*x(:,1)))+...
   sigmoid(x(:,2),-0.25).*(cos(pi/2.*x(:,2)).*cos(pi/2.*x(:,1)));

c= 2;
r= -4;
DIR1=@(x)c^2*(-r+x(:,2));
D1DIR1=@(x)zeros(size(x(:,2)));
D2DIR1=@(x)ones(size(x(:,2)));
DIR2=@(x)(-r+x(:,1));
D1DIR2=@(x)ones(size(x(:,1)));
D2DIR2=@(x)zeros(size(x(:,1)));
 
Dir_Handle=@(x,varargin)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,varargin)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];

H_Handle=@(x,flag,varargin)...
    [H1(x) H2(x)];

HEx_Handle=@(x,flag,varargin)...
    ExactSol(x,r,H_Handle);

% - curl h
VH_Handle=@(x,varargin)...
    DIR1(x)*H1(x) +DIR2(x)*H2(x);

% Load mesh from file
Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
% Mesh.BdFlags(Loc) = [-2,-2,-2,-1,-2,-1,-1,-1];
Mesh.BdFlags(Loc) = [-1,-1,-1,-2,-1,-2,-2,-2];

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
    
    NewMesh = init_LEB(NewMesh);
    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);
    NewMesh = add_DGData(NewMesh);

    % reference solution 
     H_ref = assemCochain_1f(NewMesh,HEx_Handle,gauleg(0,1,1),0);
     figure;
     plot_W1F(H_ref,NewMesh);

    
    % stabelized standard scheme
    Lie_Pen = assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LiePennW1F,Dir_Handle,gauleg(0,1,5),0.5);
    Lie_p = assemMat_W1F(NewMesh,@STIMA_ContrRot,Dir_Handle,P7O6());
    Lie_p = Lie_p + assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LieW1F,Dir_Handle,gauleg(0,1,5));
    Lie_p = Lie_p + Lie_Pen ;
    Lie_p = Lie_p - assemMat_Bnd_LieW1F(NewMesh,[-2],@STIMA_Bnd_LieW1F,Dir_Handle,gauleg(0,1,5));
    nonInFlow=find(NewMesh.BdFlags ~= -1);
    nonOutFlow=find(NewMesh.BdFlags ~= -2);
    FreeDofs=nonInFlow;

    H = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    H(FreeDofs) = 0;
    L = zeros(size(H));
    L = assemLoad_Bnd_LieW1F(NewMesh,[-1],L,gauleg(0,1,10),H_Handle, Dir_Handle,0,0);
    L = L-Lie_p*H;
    H(FreeDofs) = Lie_p(FreeDofs,FreeDofs)\L(FreeDofs);
    
    L2ErrH(j,1) = L2Err_W1F(NewMesh,H,P3O2(),HEx_Handle,0,0); 
    JPen(j,1) = sqrt((H-H_ref)'*Lie_Pen*(H-H_ref)); 
    
    plot_Norm_W1F(H-H_ref,NewMesh); colorbar;
    %     figure;
    %     plot_W1F(H,NewMesh);

    % interpolation
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());
    TopGrad = assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot = assemMat_TopRot(NewMesh);      % topological Rotation
    ContrOne = assemMat_Contr1f(NewMesh,Dir_Handle);  % contraction of one forms
    V = Dir_Handle(NewMesh.Coordinates);
    ContrTwo = assemMat_Contr2f(NewMesh,V);   % contraction of two form
    Lie_i = M*ContrTwo*TopRot+M*TopGrad*ContrOne;

    H = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    H(FreeDofs) = 0;
    cVUBd = assemBndCochain_0f(NewMesh,[-2 -3],VH_Handle);
    L = -M*TopGrad*cVUBd-Lie_i*H;
    H(FreeDofs) = (Lie_i(FreeDofs,FreeDofs))\L(FreeDofs);
    
    L2ErrH(j,2) = L2Err_W1F(NewMesh,H,P3O2(),HEx_Handle,0,0)
    JPen(j,2) = sqrt((H-H_ref)'*Lie_Pen*(H-H_ref))
    %colorbar;
    %figure; plot_W1F(H,Mesh);
    
    % exact integration
    %H_test = assemLieAdCochain_1f(NewMesh,H_Handle,Dir_Handle,JvDir_Handle,gauleg(0,1,5),0);
    %plot_Norm_W1F(H_test,NewMesh);
    %colorbar;
    %plot_Norm_W1F(H_test-H_ref,NewMesh); colorbar;
    %L2ErrH(j,3) = L2Err_W1F(NewMesh,H_test,P7O6(),HEx_Handle,0,0)
    %JPen(j,3) = sqrt((H_test-H_ref)'*Lie_Pen*(H_test-H_ref))
     
    %figure; plot_W1F(H_test,NewMesh); 
    
%     colorbar;
%     figure; plot_W1F(H_test,NewMesh);

end%  Mesh refinement

figure;
hold on;
plot(mw,L2ErrH(:,1),'rx-',...
    mw,L2ErrH(:,2),'bx-',...
    mw,L2ErrH(:,3),'gx-','Linewidth',2,'MarkerSize',8);
grid('on');
plot(mw, L2ErrH(1,1)/mw(1)*mw,'k--',mw, L2ErrH(1,1)/sqrt(mw(1))*sqrt(mw),'k--');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf L^2-error}','FontSize',12,'FontWeight','bold');
legend('Stabelized','Interpolation','TimeIntegration','h','h^{-0.5}','Location','Northwest');

slope=[diff(log(L2ErrH(:,1)))./diff(log(mw)) ...
    diff(log(L2ErrH(:,2)))./diff(log(mw)) ...
    diff(log(L2ErrH(:,3)))./diff(log(mw))]

figure;
hold on;
plot(mw,JPen(:,1),'rx-',...
    mw,JPen(:,2),'bx-',...
    mw,JPen(:,3),'gx-','Linewidth',2,'MarkerSize',8);
grid('on');
plot(mw, JPen(1,1)/mw(1)*mw,'k--',mw, JPen(1,1)/sqrt(mw(1))*sqrt(mw),'k--');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf JumpPen-error}','FontSize',12,'FontWeight','bold');
legend('Stabelized','Interpolation','TimeIntegration','h','h^{-0.5}','Location','Northwest');

slope=[diff(log(JPen(:,1)))./diff(log(mw)) ...
    diff(log(JPen(:,2)))./diff(log(mw)) ...
    diff(log(JPen(:,3)))./diff(log(mw))]

figure;
hold on;
plot(mw,L2ErrH(:,1),'rx-',...
    mw,L2ErrH(:,2),'bx-','Linewidth',2,'MarkerSize',8);
grid('on');
plot(mw, L2ErrH(1,1)/mw(1)*mw,'k--',mw, L2ErrH(1,1)/sqrt(mw(1))*sqrt(mw),'k--');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf L^2-error}','FontSize',12,'FontWeight','bold');
legend('Stabelized','Interpolation','h','h^{-0.5}','Location','Northwest');



function h = ExactSol(x,r,Hhandle)
n = size(x,1);
h = zeros(n,2);
if r == -2;
    for i =1:n
        d1 = sqrt(-8+4*x(i,1)+x(i,1)^2-16*x(i,2)-4*x(i,2)^2);
        d2 = sqrt(13-4*x(i,1)-x(i,1)^2+16*x(i,2)+4*x(i,2)^2);
        if abs(x(i,1)-2*(1+x(i,2)))<eps
            xloc = [0,-1];
            T = 1/(2+x(i,2));
        else
            if isreal(d1) && d1>=1
                xloc = [-2+d1,-1];
                T = (-2-2*x(i,2)+x(i,1))/xloc(1);
            else
                if isreal(d2) && d2>=2
                    xloc = [-1,-2+d2/2];
                    T = (2+2*x(i,2)-x(i,1))/(3+2*xloc(2));
                end
            end
        end
        DX = 1/2*[T+1/T 2*(T-1/T);  1/2*(T-1/T) T+1/T];
        h(i,:) = Hhandle(xloc)*DX;
    end
end
if r ==-4
    for i =1:n
        d1 = sqrt(-12+8*x(i,1)+x(i,1)^2-32*x(i,2)-4*x(i,2)^2);
        d2 = sqrt(57-8*x(i,1)-x(i,1)^2+32*x(i,2)+4*x(i,2)^2);
        if abs(x(i,1)-2*(2+x(i,2)))<eps
            xloc = [2,-1];
            T = 3/(4+x(i,2));
        else
            if isreal(d1) && d1>=3
                xloc = [-4+d1,-1];
                T = (-4-2*x(i,2)+x(i,1))/(-2+xloc(1));
            else
                if isreal(d2) && d2>=6
                    xloc = [-1,-4+d2/2];
                    T = (4 +2*x(i,2)-x(i,1))/(5+2*xloc(2));
                end
            end
        end
        DX = 1/2*[T+1/T 2*(T-1/T);  1/2*(T-1/T) T+1/T];
        h(i,:) = Hhandle(xloc)*DX;
    end
end

function y=sigmoid(x,x0)
n=max(size(x));
y = zeros(size(x));
for i = 1:n
    if x(i)-x0<=0.5 &&  x(i)-x0>=-0.5
        %y(i)=1/2*(1+ 10^6*(35/1600000*(x(i)-x0)-35/16000*(x(i)-x0).^3+21/160*(x(i)-x0).^5-50/16*(x(i)-x0).^7));
        y(i)=(0.5+ 2.1875*(x(i)-x0)-8.75*(x(i)-x0).^3+21*(x(i)-x0).^5-20*(x(i)-x0).^7);
    else if x(i)-x0>0.5
            y(i)=1;
        end
    end
end
