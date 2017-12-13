%main_ContourExample

NREFS = 5;              %Number of red refinement steps
F_HANDLE = @f_LShap;    %Right hand side of source term
GD_HANDLE = @g_D_LShap; %Dirichlet boundary data
GN_HANDLE = @g_N_LShap; %Neumann boundary data

% Initializing mesh

Mesh = load_Mesh('Coord_LShap.dat', 'Elem_LShap.dat'); %Sets up standard L-mesh
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
Mesh = add_Edges(Mesh); %For each elements you get all edges
Loc = get_BdEdges(Mesh); %Extract the Borderedges of the mesh
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); %Initialize the Flags on the border
Mesh.BdFlags(Loc) = [-1 -2 -7 -7 -3 -4];
for i = 1:NREFS,
    Mesh = refine_REG(Mesh); %Refining the mesh
end
Mesh = add_Edge2Elem(Mesh); %Defining which elements is on each side of the edge

%Assemble stiffness matrix and load vector

A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE); %Making the stiffnesmatrix
L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE); %Making the righthand side of eqn

%Incorporate Neumann boundary data

L=assemNeu_LFE(Mesh,-1:-1:-4,L,gauleg(0,1,4),GN_HANDLE);

% Incorporate Dirichlet boundary data

[U,DF] = assemDir_LFE(Mesh,-7,GD_HANDLE);
L = L - A*U; %Put the already known U's on the right side

%Solves the system for the real unknowns

U(DF) = A(DF,DF)\L(DF);

%Contourplot of the solution

contour_LFE(U,Mesh);