% produces plots from the data obtained by running main_exp1.m

%   Copyright 2008-2008 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

data = load('experiment1.mat');
a = data.exp1_diffrange;

d=getData(2);


%%%%%%%%%%%%%%%%
% Mesh plots and solution plots
%%%%%%%%%%%%%%%%

% produce a plot of the exact solution

Mesh = data.exp1_meshes(end);
sol = d.U_EX_Handle(Mesh.Coordinates,1,1e-10);
plot_LFE(sol,Mesh);colorbar;
FileName = 'Exact_sol.eps';
print('-depsc2', FileName); 

% produce a mesh plot of a coarse mesh
Mesh = data.exp1_meshes(1);
plot_Mesh(Mesh,'as');
FileName = 'Exp1Mesh.eps';
print('-depsc2', FileName); 

% plot the solution for the different methods and different diffusion
% constants
Mesh = data.exp1_meshes(end);
U = data.exp1_solStGal(end,2).sol;
U_up = data.exp1_solUP(end,2).sol;
U_supg = data.exp1_solSUPG(end,2).sol;
U_fvol = data.exp1_solFVol(end,2).sol;

% plot_LFE(U,Mesh);colorbar;
% plot_LFE(U_up,Mesh);colorbar;
% plot_LFE(U_supg,Mesh);colorbar;
plotLine_LFE(U,Mesh,[0 0],[1 0]);
FileName = ['galerkin' num2str(a(2)) '.eps'];
print('-depsc2', FileName); 
plotLine_LFE(U_up,Mesh,[0 0],[1 0]);
FileName = ['upwind' num2str(a(2)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_supg,Mesh,[0 0],[1 0]);
FileName = ['supg' num2str(a(2)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_fvol,Mesh,[0 0],[1 0]);
FileName = ['fvol' num2str(a(2)) '.eps'];
print('-depsc2', FileName);


U = data.exp1_solStGal(end,4).sol;
U_up = data.exp1_solUP(end,4).sol;
U_supg = data.exp1_solSUPG(end,4).sol;
U_fvol = data.exp1_solFVol(end,4).sol;

% plot_LFE(U,Mesh);colorbar;
% plot_LFE(U_up,Mesh);colorbar;
% plot_LFE(U_supg,Mesh);colorbar;
plotLine_LFE(U,Mesh,[0 0],[1 0]);
FileName = ['galerkin' num2str(a(4)) '.eps'];
print('-depsc2', FileName); 
plotLine_LFE(U_up,Mesh,[0 0],[1 0]);
FileName = ['upwind' num2str(a(4)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_supg,Mesh,[0 0],[1 0]);
FileName = ['supg' num2str(a(4)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_fvol,Mesh,[0 0],[1 0]);
FileName = ['fvol' num2str(a(4)) '.eps'];
print('-depsc2', FileName);

U = data.exp1_solStGal(end,6).sol;
U_up = data.exp1_solUP(end,6).sol;
U_supg = data.exp1_solSUPG(end,6).sol;
U_fvol = data.exp1_solFVol(end,6).sol;

plot_LFE(U,Mesh);colorbar;
FileName = ['fullgalerkin' num2str(a(6)) '.eps'];
print('-depsc2', FileName); 
% plot_LFE(U_up,Mesh);colorbar;
% plot_LFE(U_supg,Mesh);colorbar;
plotLine_LFE(U,Mesh,[0 0],[1 0]);
FileName = ['galerkin' num2str(a(6)) '.eps'];
print('-depsc2', FileName); 
plotLine_LFE(U_up,Mesh,[0 0],[1 0]);
FileName = ['upwind' num2str(a(6)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_supg,Mesh,[0 0],[1 0]);
FileName = ['supg' num2str(a(6)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_fvol,Mesh,[0 0],[1 0]);
FileName = ['fvol' num2str(a(6)) '.eps'];
print('-depsc2', FileName);

%%%%%%%%%%%%%%%%%%%%
% plots of the errors
%%%%%%%%%%%%%%%%%%%%

h = data.exp1_h;
NREFS = length(h);

% plot the L^2 error of the 3 methods for different diffusion constants
clear err;
err(:,1) = data.exp1_errStGal(end,:);
err(:,2) = data.exp1_errUPQuad(end,:);
err(:,3) = data.exp1_errSUPG(end,:);
err(:,4) = data.exp1_errFVol(end,:);

fig = figure('Name','Discretization error');
plot(a,err(:,1),'go-',a,err(:,2),'bo-',a,err(:,3),'ro-',a,err(:,4),'co-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf Diffusion Constant}');
ylabel('{\bf Error}');
legend('L^2-error u','L^2-error u (up)','L^2-error u (supg)','L^2-error u (fvol)','Location','NorthEast');

FileName = 'err_diff.eps';
print('-depsc2', FileName);


% plot the L^2 error of the 3 methods for different refinement steps
clear err;
err(:,1) = data.exp1_errStGal(:,1);
err(:,2) = data.exp1_errUPQuad(:,1);
err(:,3) = data.exp1_errSUPG(:,1);
err(:,4) = data.exp1_errFVol(:,1);

fig = figure('Name',['Discretization error (a = ',num2str(a(1)),')']);
plot(h,err(:,1),'go-',h,err(:,2),'bo-',h,err(:,3),'ro-',h,err(:,4),'co-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf Error}');
legend('L^2-error u','L^2-error u (up)','L^2-error u (supg)','L^2-error u (fvol)','Location','NorthEast');
p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,1)),1);
add_Slope(gca,'Northwest',p(1),'g-');
p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,2)),1);
add_Slope(gca,'East',p(1),'b-');
p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,3)),1);
add_Slope(gca,'Southeast',p(1),'r-');
p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,4)),1);
add_Slope(gca,'North',p(1),'c-');

FileName = ['discr_err' num2str(a(1)) '.eps'];
print('-depsc2', FileName);


err(:,1) = data.exp1_errStGal(:,end);
err(:,2) = data.exp1_errUPQuad(:,end);
err(:,3) = data.exp1_errSUPG(:,end);
err(:,4) = data.exp1_errFVol(:,end);

fig = figure('Name',['Discretization error (a = ',num2str(a(end)),')']);
plot(h,err(:,1),'go-',h,err(:,2),'bo-',h,err(:,3),'ro-',h,err(:,4),'co-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf Error}');
legend('L^2-error u','L^2-error u (up)','L^2-error u (supg)','L^2-error u (fvol)','Location','NorthEast');
p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,1)),1);
add_Slope(gca,'West',p(1),'g-');
p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,2)),1);
add_Slope(gca,'East',p(1),'b-');
p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,3)),1);
add_Slope(gca,'Southwest',p(1),'r-');
p = polyfit(log(h(1:NREFS)),log(err(1:NREFS,4)),1);
add_Slope(gca,'North',p(1),'c-');

FileName = ['discr_err' num2str(a(end)) '.eps'];
print('-depsc2', FileName);

