% LehrFEM MATLAB script for computing Dirichlet eigenvalues of Laplacian
% on a unit disc domain.

GD_HANDLE = @(x,varargin)zeros(size(x,1),1); % Zero Dirichlet data
H0 =[ .25 .2 .1 .05 .02 .01 0.005]'; % target mesh widths
NRef = length(H0); % Number of refinement steps

% Variables for mesh widths and eigenvalues
M_W = zeros(NRef,1); lmax = M_W; lmin = M_W;

% Main refinement loop
for iter = 1:NRef 
  
% Set parameters for mesh
  C = [0 0];             % Center of circle
  R = 1;                 % Radius of circle 
  BBOX = [-1 -1; 1 1];   % Bounding box
  DHANDLE = @dist_circ;  % Signed distance function
  HHANDLE = @h_uniform;  % Element size function
  FIXEDPOS = [];         % Fixed boundary vertices of the mesh
  DISP = 0;              % Display flag
  
% Mesh generation
  Mesh = init_Mesh(BBOX,H0(iter),DHANDLE,HHANDLE,FIXEDPOS,DISP,C,R);  
  Mesh = add_Edges(Mesh);  % Provide edge information
  Loc = get_BdEdges(Mesh); % Obtain indices of edges on $\partial\Omega$
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;  % Flag boundary edges 
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  M_W(iter) = get_MeshWidth(Mesh); % Get mesh width
  
  fprintf('Mesh on level %i: %i elements, h = %f\n',iter,size(Mesh,1),M_W(iter));
% Assemble stiffness matrix and mass matrix
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE,P7O6());
  M = assemMat_LFE(Mesh,@MASS_LFE,P7O6());
% Incorporate Dirichlet boundary data (nothing to do here)
  [U,FreeNodes] = assemDir_LFE(Mesh,-1,GD_HANDLE);
  A = A(FreeNodes,FreeNodes);
  M = M(FreeNodes,FreeNodes);
  
% Use MATLAB's built-in \texttt{eigs}-function to compute the
% extremal eigenvalues, see \ncseref[Sect.]{sec:KrylovEigensolvers}.
  NEigen = 6;
  d = eigs(A,M,NEigen,'sm'); lmin(iter) = min(d);
  d = eigs(A,M,NEigen,'lm'); lmax(iter) = max(d);
end

figure; plot(M_W,lmin,'b-+',M_W,lmax,'r-*'); grid on;
set(gca,'XScale','log','YScale','log','XDir','reverse');
title('\bf Eigenvalues of Laplacian on unit disc');
xlabel('{\bf mesh width h}','fontsize',14);
ylabel('{\bf generalized eigenvalues}','fontsize',14);
legend('\lambda_{min}','\lambda_{max}','Location','NorthWest')
p = polyfit(log(M_W),log(lmax),1);
add_Slope(gca,'east',p(1));

print -depsc2 '../../../Slides/NPDEPics/geneigdisklfe.eps';

