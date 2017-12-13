function main_SLPureCircErrorTestW1F()
%  Semi-Lagrange- version for MHD

%   Copyright 2011 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

addpath('./Lib')
clear Mesh;
%clear all;
NREFS_init =0;                             % Number of initial uniform red refinements
NREFS = 2;                                   % Number of uniform red refinements
CFLs = [0.4]  ;                              %  CFL-numbers 0.5,1,1.5,2,2.5,
T0 = 0;                                         % start time 
T1= pi;                                         % stop time
locIntegrate = 'Implicit Euler';   % local integrator: 'Midpoint Rule', 'Explict Euler', 'Implicit Euler'
nsteps_loc =1                             % local timesteps

% allocate memory
tau =zeros(NREFS,size(CFLs,2));
steps =zeros(NREFS,size(CFLs,2));
error =zeros(NREFS,size(CFLs,2));
h =zeros(NREFS,1);

% inital data, Gaussian times something
lambda=0.1;
H1=@(x)exp(-1/(2*lambda^2)*((x(:,1)-0.5).^2+(x(:,2)).^2));
H2=@(x)exp(-1/(2*lambda^2)*((x(:,1)-0.5).^2+(x(:,2)).^2));
H_Handle=@(x,flag,t,varargin)[H1(x) H2(x)];

% zero right hand side
F_Handle=@(x,flag,t,varargin)zeros(size(x));

% mesh generation
C = [0 0];                             % Center of circle
R = 1;                                   % Radius of circle
BBOX = [-1 -1; 1 1];           % Bounding box
H0 = 0.5;                             % Initial mesh width
DHANDLE = @dist_circ;      % Signed distance function
HHANDLE = @h_uniform;  % Element size function
FIXEDPOS = [];                    % Fixed boundary vertices of the mesh
DISP = 1;                              % Display flag

Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP,C,R);
Mesh = add_Edges(Mesh);
%plot_Mesh(Mesh,'as');

% velocity field, set to 0 at boundary
v1=1;
v2=1;
DIR1=@(x)x(:,2);
D1DIR1=@(x)zeros(size(x,1),1);
D2DIR1=@(x)ones(size(x,1),1);
DIR2=@(x)-x(:,1);
D1DIR2=@(x)-ones(size(x,1),1);
D2DIR2=@(x)zeros(size(x,1),1);

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];

% Add edge data structure and boundary flags
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

%  refine mesh
for i=1:NREFS_init
      Mesh = refine_REG(Mesh,DHANDLE,C,R);
      %plot_Mesh(Mesh,'a');
end

NewMesh = Mesh;

for j=1:NREFS
      
    NewMesh = refine_REG(NewMesh,DHANDLE,C,R);
    
    % add edge2elements and patches data
    NewMesh = init_LEB(NewMesh);
    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);

    NewMesh = add_DGData(NewMesh);

    
    % mesh width
    h(j) = get_MeshWidth(NewMesh);
    for k = 1:size(CFLs,2)
        
        % mesh width, time step size and CFL-number
        tau(j,k) = CFLs(k)*h(j)^1/norm([v1,v2]);
        steps(j,k) = ceil((T1-T0)/tau(j,k));
        tau(j,k) = tau(j,k)-(T0+steps(j,k)*tau(j,k)-T1)/steps(j,k);
        
        
        % discrete pullback via interpolation
        % pbVm = trace_vertices_WodeSol(NewMesh,Dir_Handle,-tau);
        pbVm = trace_vertices_RotationMidpoint(NewMesh,-tau(j,k),nsteps_loc, locIntegrate);
        pbVp = trace_vertices_RotationMidpoint(NewMesh,tau(j,k),nsteps_loc, locIntegrate);
        %pbVm = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,-h);
        
        %plot_RotMesh(NewMesh,pbVm,'a');
        P_i = assemMat_SemiLag_W1Fcheat(NewMesh, pbVm);
        
        defMesh = NewMesh;
        defMesh.Coordinates = pbVp(:,[1 2]);
        plot_transportMesh(NewMesh,pbVp,'pe ');
        intersec = aff_elems4(NewMesh,pbVp);
        plot_MeshPatch(NewMesh,defMesh,pbVp,intersec,' ');
        P_p = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);

        % Assemble MASS matrix
        %M = assemMat_W1F(NewMesh,@MASS_W1F,@(x,varargin)1, P3O3());
        M = assemMat_MASSQuad_W1F_patches(NewMesh, defMesh,intersec);
       
        % inital data
        H_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
        H_oldp=[H_old];
        plot_Norm_W1F(H_old(:,1),NewMesh);
        colorbar; CLim1=get(gca,'CLim');
        plot_Norm_W1F(H_oldp(:,1),NewMesh);
        colorbar; CLim1=get(gca,'CLim');
        %figure; plot_W1F(H_old(:,1),NewMesh);

        % timestepping
        time(1) = 0;
        for i = 1:steps(j,k)
            [i steps(j,k)]
            time(i+1) = tau(j,k)+time(i);

            %Dirichlet data
            %H_new=zeros(size(H_old));

            %righthand side
            %L = assemCochain_1f(NewMesh, F_Handle,gauleg(0,1,10),i*tau(j,k));

            %iterate system
            H_old = (P_i*H_old);
            H_oldp = M\(P_p*H_oldp);
            
            plot_Norm_W1F(H_old(:,1),NewMesh); 
            set(gca,'CLim',CLim1); colorbar;
            plot_Norm_W1F(H_oldp(:,1),NewMesh); 
            set(gca,'CLim',CLim1); colorbar;

            %Update vectors
            %H_old=H_new;

        end % timestep

        plot_Norm_W1F(H_old(:,1),NewMesh);
        %rectangle('Position',[-0.5,-0.5,1,1],'Curvature',[1,1])
        set(gca,'CLim',CLim1); colorbar;
       
        error(j,k) = L2Err_W1F(NewMesh,H_old,P7O6(),H_Handle,0,0);
    end %CFLS
end % mesh

data = [h error]