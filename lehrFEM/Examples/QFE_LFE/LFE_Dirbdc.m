% Initialize constants anf functions
F_HANDLE = @f_LShap;     % Right hand side source term
GD_HANDLE = @g_D_LShap;  % Dirichlet boundary data
% Load mesh
Mesh = load_Mesh('meshvert.dat','meshel.dat'); 
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
Mesh = add_Edges(Mesh);                                
Loc = get_BdEdges(Mesh); % Obtain indices of edges on \Blue{$\partial\Omega$}
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = 1;   % Set flag '1' for edges on the boundary
  
% Assemble Galerkin (stiffness) matrix and right hand side (load) vector for linear Lagrangian FE
A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
% Note: the arguments \texttt{0,1,2} are passed on to \texttt{F\_HANDLE} through
% \texttt{varargin} and they are specific to this particular example.
phi = assemLoad_LFE(Mesh,P7O6(),F_HANDLE,0,1,2);

% Incorporate Dirichlet boundary data for vertices adjacent to edges
% carrying flag $1$. \texttt{U} contains nodal values for Dirichlet
% boundary data, \texttt{FreeDofs} contains numbers of interior nodes.
% See Ex.~\ref{ex:dirbdcimp} for further explanations.
[U,FreeDofs] = assemDir_LFE(Mesh,[1],GD_HANDLE);
phi = phi - A*U; % \label{db:1}
  
% Solve the linear system with matrix \Blue{$\VA_0$} from \eqref{eq:fullLSE}
U(FreeDofs) = A(FreeDofs,FreeDofs)\phi(FreeDofs);

% Plot solution
plot_LFE(U,Mesh);colorbar;
plotLine_LFE(U,Mesh,[0 0] ,[1 1]);
