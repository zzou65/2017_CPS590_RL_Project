clear all
close all
cd ./lehrFEM/LehrFEM/;
startup;
cd ../../

AD2D_Global_Vars;

n_mesh_refine = 5;
nDim = 2;
data.boundtype =[0 0 0 -1]; %b, l, r, t
data.BC_Handle = @(x,varargin) zeros(size(x, 1), 1);
data.F_Handle  = @(x, elemflag, xc, yc, sigma) ...
                    exp(-((x(:,1)-xc).^2+(x(:,2)-yc).^2)./(0.5*sigma^2));
data.V_Handle_L = @(x,varargin) AD2D_V_Handle(x, 1, 0.5); 
data.V_Handle_R = @(x,varargin) AD2D_V_Handle(x, 2, 0.5); 
data.Coordinates=[0 0; 1 0; 1 1; 0 1];
data.Elements=[1 2 4;2 3 4];
DataObj = data;

% prepare Mesh obj
Mesh.Coordinates =data.Coordinates; 
Mesh.Elements = data.Elements;
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = data.boundtype;
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
for i=1:n_mesh_refine
    Mesh = refine_REG(Mesh);
end
Mesh = add_Edge2Elem(Mesh);
MeshObj = Mesh;
EleCenterCoords = computeEleCenters(MeshObj);
NEle = size(MeshObj.Elements,1);
NNode = size(MeshObj.Coordinates, 1);

%%% DBC 
[U_s,FreeDofs] = assemDir_LFE(MeshObj,[-1],data.BC_Handle);
U_w_BC = U_s;

%%% Matrices
M = assemMat_LFE(MeshObj,@MASS_LFE,P7O6());
% test M matrix:
vec = ones(NNode, 1);
vec_integral = vec'*M*vec;
if(abs(vec_integral - 1) < 1e-3)
    fprintf('M check passed ...\n')
end

S_Diff_Const = 2.5e-2;
S_Diff_Field = S_Diff_Const * ones(NEle, 1);
DS = AD2D_assemMat_LFE(MeshObj,@STIMA_Lapl_LFE, S_Diff_Field, P3O3());
C_Diff_Const = 5.0e-2;
C_Diff_Field = C_Diff_Const * ones(NEle, 1);
DC = AD2D_assemMat_LFE(MeshObj,@STIMA_Lapl_LFE, C_Diff_Field, P3O3());
AL = assemMat_LFE(MeshObj, @STIMA_Conv_LFE, DataObj.V_Handle_L, P7O4());
AR = assemMat_LFE(MeshObj, @STIMA_Conv_LFE, DataObj.V_Handle_R, P7O4());

%%% setup sources
stepSize = 1/8;
FCenters = [1/4, 1/4; 3/4, 1/4];
FSigmas = [1/10; 1/10];
xgv = 0:stepSize:1;
ygv = xgv;
[X,Y] = meshgrid(xgv,ygv);
X=X(:);
Y=Y(:);
add_centers = [X, Y];
n_add = length(X);
add_sigmas = 1/10 * ones(n_add, 1);
FCenters = [FCenters; add_centers];
FSigmas = [FSigmas; add_sigmas];
FVecs = zeros(NNode, length(FSigmas));
for i=1:length(FSigmas)
    if(mod(i, 10) == 1)
        fprintf('Makeing the %d-th RHS...\n', i)
    end
    F = assemLoad_LFE(MeshObj,P3O3(), data.F_Handle, FCenters(i, 1), FCenters(i, 2), FSigmas(i));
    small_idx = find(abs(F) < 1e-12);
    F(small_idx) = 0.0;
    FVecs(:, i) = F;
end
FC = -FVecs(:, 1:2);
FS = FVecs(:, 3:end);
SPos = add_centers;

%%% reward
cost_thre_low = -0.05; 
cost_thre_hig =  0.02;
costL = 2e-2;
costR = 2e-2;

%%% algorithmic parameters
theta = 0.5;
n_sub_steps = 6;
dt = 1 / n_sub_steps;

%%% Reset system state
AD2D_ResetState();

%%% For image representation
N_Pixels_Per_Dim = 20;
Image_Projection_Mat = AD2D_MakeImgProjMat( N_Pixels_Per_Dim );

N_Add_Feature = 3;
N_Pixels_Per_Dim_for_Feature = 4;
Image_Projection_Mat_for_Feature = AD2D_MakeImgProjMat( N_Pixels_Per_Dim_for_Feature );
If_Use_AddFeature = 1;
if(If_Use_AddFeature == 1)
    N_Feature = N_Add_Feature + N_Pixels_Per_Dim_for_Feature ^ 2;
else
    N_Feature = N_Pixels_Per_Dim_for_Feature ^ 2;
end
    
% reward_low = -0.12;
% reward_hig =  0.10;
reward_low = -0.12;
reward_hig =  0.12;
state_low = -0.1;
state_hig = 0.15;

% %%% plot F field:w
 
% figure(1)
% r = randi([1 n_add],1,1);
% subplot(1,2,1)
% F_field = 0;
% F_field = F_field - data.F_Handle(EleCenterCoords,[], ...
%         FCenters(1, 1), FCenters(1, 2), FSigmas(1));
% F_field = F_field - data.F_Handle(EleCenterCoords,[], ...
%         FCenters(2, 1), FCenters(2, 2), FSigmas(1));
% F_field = F_field + data.F_Handle(EleCenterCoords,[], ...
%         FCenters(r + 2, 1), FCenters(r + 2, 2), FSigmas(r + 2));    
% plotSurf(EleCenterCoords, F_field, [-1, 1], 1)
% subplot(1,2,2)
% F_field = 0;
% F_field = F_field + FC(:, 1);
% F_field = F_field + FC(:, 2);
% F_field = F_field + FS(:, r);   
% plotSurf(Mesh.Coordinates, F_field, [], 1)
