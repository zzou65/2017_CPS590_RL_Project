% produces plots from the data obtained by running main_exp2.m

%   Copyright 2008-2008 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

data = load('experiment2.mat');
a = data.exp2_diffrange;

%%%%%%%%%%%%%%%%
% Mesh plots and solution plots
%%%%%%%%%%%%%%%%

% produce a mesh plot of a coarse mesh
Mesh = data.exp2_meshes(1);
plot_Mesh(Mesh,'as');
FileName = 'Exp2Mesh.eps';
print('-depsc2', FileName); 

% plot the solution for the different methods and different diffusion
% constants
Mesh = data.exp2_meshes(end);

%%%%%%
U = data.exp2_solStGal(end,1).sol;
U_up = data.exp2_solUP(end,1).sol;
U_supg = data.exp2_solSUPG(end,1).sol;
U_fvol = data.exp2_solFVol(end,1).sol;

% full region plots
plot_LFE(U,Mesh);colorbar;
FileName = ['fullgal_exp2' num2str(a(1)) '.eps'];
print('-depsc2', FileName); 
plot_LFE(U_up,Mesh);colorbar;
FileName = ['fullupw_exp2' num2str(a(1)) '.eps'];
print('-depsc2', FileName); 
plot_LFE(U_supg,Mesh);colorbar;
FileName = ['fullsupg_exp2' num2str(a(1)) '.eps'];
print('-depsc2', FileName); 
plot_LFE(U_fvol,Mesh);colorbar;
FileName = ['fullfvol_exp2' num2str(a(1)) '.eps'];
print('-depsc2', FileName); 

% line plots
plotLine_LFE(U,Mesh,[0 1],[1 0]);
FileName = ['galerkin_exp2' num2str(a(2)) '.eps'];
print('-depsc2', FileName); 
plotLine_LFE(U_up,Mesh,[0 1],[1 0]);
FileName = ['upwind_exp2' num2str(a(2)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_supg,Mesh,[0 1],[1 0]);
FileName = ['supg_exp2' num2str(a(2)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_fvol,Mesh,[0 1],[1 0]);
FileName = ['fvol_exp2' num2str(a(2)) '.eps'];
print('-depsc2', FileName);

%%%%%%
U = data.exp2_solStGal(end,4).sol;
U_up = data.exp2_solUP(end,4).sol;
U_supg = data.exp2_solSUPG(end,4).sol;
U_fvol = data.exp2_solFVol(end,4).sol;

% full region plots
plot_LFE(U,Mesh);colorbar;
FileName = ['fullgal_exp2' num2str(a(4)) '.eps'];
print('-depsc2', FileName); 
plot_LFE(U_up,Mesh);colorbar;
FileName = ['fullupw_exp2' num2str(a(4)) '.eps'];
print('-depsc2', FileName); 
plot_LFE(U_supg,Mesh);colorbar;
FileName = ['fullsupg_exp2' num2str(a(4)) '.eps'];
print('-depsc2', FileName); 
plot_LFE(U_fvol,Mesh);colorbar;
FileName = ['fullfvol_exp2' num2str(a(4)) '.eps'];
print('-depsc2', FileName); 

%line plots
plotLine_LFE(U,Mesh,[0 1],[1 0]);
FileName = ['galerkin_exp2' num2str(a(4)) '.eps'];
print('-depsc2', FileName); 
plotLine_LFE(U_up,Mesh,[0 1],[1 0]);
FileName = ['upwind_exp2' num2str(a(4)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_supg,Mesh,[0 1],[1 0]);
FileName = ['supg_exp2' num2str(a(4)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_fvol,Mesh,[0 1],[1 0]);
FileName = ['fvol_exp2' num2str(a(4)) '.eps'];
print('-depsc2', FileName);

%%%%%%
U = data.exp2_solStGal(end,6).sol;
U_up = data.exp2_solUP(end,6).sol;
U_supg = data.exp2_solSUPG(end,6).sol;
U_fvol = data.exp2_solFVol(end,6).sol;

% full region plots
plot_LFE(U,Mesh);colorbar;
FileName = ['fullgal_exp2' num2str(a(6)) '.eps'];
print('-depsc2', FileName); 
plot_LFE(U_up,Mesh);colorbar;
FileName = ['fullupw_exp2' num2str(a(6)) '.eps'];
print('-depsc2', FileName); 
plot_LFE(U_supg,Mesh);colorbar;
FileName = ['fullsupg_exp2' num2str(a(6)) '.eps'];
print('-depsc2', FileName); 
plot_LFE(U_fvol,Mesh);colorbar;
FileName = ['fullfvol_exp2' num2str(a(6)) '.eps'];
print('-depsc2', FileName); 

%line plots
plotLine_LFE(U,Mesh,[0 1],[1 0]);
FileName = ['galerkin_exp2' num2str(a(6)) '.eps'];
print('-depsc2', FileName); 
plotLine_LFE(U_up,Mesh,[0 1],[1 0]);
FileName = ['upwind_exp2' num2str(a(6)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_supg,Mesh,[0 1],[1 0]);
FileName = ['supg_exp2' num2str(a(6)) '.eps'];
print('-depsc2', FileName);
plotLine_LFE(U_fvol,Mesh,[0 1],[1 0]);
FileName = ['fvol_exp2' num2str(a(6)) '.eps'];
print('-depsc2', FileName);


% clear memory

clear all;

