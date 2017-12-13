function main_SLPureCircErrorRate_hp(NREFS, CFLs, locIntegrate, Problem,timestepchoice,polydegree)
%  Semi-Lagrange- version for linear advection with continous LFE

%   Copyright 2011 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

addpath('./Lib')
clear Mesh;
%clear all;

% parameters
NREFS_init = 0;                            % Number of initial uniform red refinements
%NREFS = 3;                                 % Number of uniform red refinements
%CFLs = [0.4]  ;                            % CFL-numbers 0.5,1,1.5,2,2.5,
T0 = 0;                                          % start time 
T1= pi/2;                                      % stop time
%locIntegrate = 'Midpoint Rule'; % local integrator: 'Midpoint Rule', 
                                                       %                  'Explict Euler',
                                                       %                  'Implicit Euler',
                                                       %                  'Trapez Rule',
                                                       %                  'Explicit Midpoint Rule'.
                                                       %                  'exact'
nsteps_loc = 1                              % local timesteps
% polydegree =2                             % polynomial degree
% Problem = 'Gauss'                      % Problemtype 'Gauss','CosQ'
% timestepchoice = -1                  % use CFL condition if timestepchoise < 0 ,
                                                        % and fixed timestep timestepchoice otherwise    
                                                        
% If we use inexact quadrature on the polygons for the assembling of the pullback
% we use the same assembling (inexact quadrature on polygons) for the mass matrix,
% otherwise we seem to have stability problems.
quadflag = 'nobary';                         % 'bary': simple barycenter quadrature; otherwise Quadrule in triangulation
QuadRule = P7O6();                     % quadrature rule on the 

% allocate memory
tau =zeros(NREFS,size(CFLs,2));
steps =zeros(NREFS,size(CFLs,2));
error =zeros(NREFS,size(CFLs,2));
h =zeros(NREFS,1);

% inital data, Gaussian times something

if strcmp(Problem,'Gauss')
    lambda=0.2;
    x0 = 0;
    y0 = 0.5;
    u=@(x)exp(-1/(2*lambda^2)*((x(:,1)-x0).^2+(x(:,2)-y0).^2));
    u_Handle=@(x,flag,t,varargin)u(x);
end

if strcmp(Problem,'CosQ')
    x0 = 0;
    y0 = 0.25;
    a = 2;
    us=@(x) pi/2*a.*sin(pi*a*sqrt((x(:,1)-x0).^2+(x(:,2)-y0).^2))./sqrt((x(:,1)-x0).^2+(x(:,2)-y0).^2);
    u=@(x) us(x).*Const_pulse_2D(x,1/a,[x0,y0]).*(x0-x(:,1));
    u_Handle=@(x,flag,t,varargin)u(x);
end

if strcmp(Problem,'CosQQ')
    x0 = 0;
    y0 = 0.25;
    a = 2;
    us=@(x) pi/4*a.*(2*sin(pi*a*sqrt((x(:,1)-x0).^2+(x(:,2)-y0).^2))+...
        sin(2*pi*a*sqrt((x(:,1)-x0).^2+(x(:,2)-y0).^2)))./sqrt((x(:,1)-x0).^2+(x(:,2)-y0).^2);
    u=@(x) us(x).*Const_pulse_2D(x,1/a,[x0,y0]).*(x0-x(:,1));
    u_Handle=@(x,flag,t,varargin)u(x);
end

sol_Handle=@(x,t,varargin)u([cos(t)*x(:,1)-sin(t)*x(:,2) sin(t)*x(:,1)+cos(t)*x(:,2)]);

% zero right hand side
F_Handle=@(x,flag,t,varargin)zeros(size(x));

% mesh generation
C = [0 0];                             % Center of circle
R = 1;                                 % Radius of circle
BBOX = [-1 -1; 1 1];                   % Bounding box
H0 = 0.5;                              % Initial mesh width
DHANDLE = @dist_circ;                  % Signed distance function
HHANDLE = @h_uniform;                  % Element size function
FIXEDPOS = [];                         % Fixed boundary vertices of the mesh
DISP = 1;                              % Display flag

Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP,C,R);
Mesh = add_Edges(Mesh);
%plot_Mesh(Mesh,'as');

% Add edge data structure and boundary flags
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1; 
Mesh.ElemFlag =zeros(size(Mesh.Elements,1),1); 

%  refine mesh
for i=1:NREFS_init
      Mesh = refine_REG(Mesh,DHANDLE,C,R);
      %plot_Mesh(Mesh,'a');
end

NewMesh = Mesh;

% Build shape functions and quadrature rules
QuadRule_1D = gauleg(0,1,2*polydegree);
Shap_1D = shap_hp([QuadRule_1D.x zeros(size(QuadRule_1D.x))],polydegree);
QuadRule_2D = Duffy(TProd(QuadRule_1D));
Shap_2D = shap_hp(QuadRule_2D.x,polydegree);

for j=1:NREFS
      
    NewMesh = refine_REG(NewMesh,DHANDLE,C,R);

    % add edge2elements, patches data
    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);

    % Assign polynomial degrees and build dof maps
    EDofs = (polydegree-1)*ones(size(NewMesh.Edges,1),1);
    if(p > 2)
        CDofs = (polydegree-1)*(polydegree-2)/2*ones(size(NewMesh.Elements,1),1);
    else
        CDofs = zeros(size(NewMesh.Elements,1),1);
    end
    Elem2Dof = build_DofMaps(NewMesh,EDofs,CDofs);
    
    % mesh width
    h(j) = get_MeshWidth(NewMesh);
    for k = 1:size(CFLs,2)
        
        % mesh width, time step size and CFL-number
        if  timestepchoice < 0
            tau(j,k) = CFLs(k)*h(j)^1/sqrt(2);
        else
            tau(j,k) = timestepchoice;
        end
        steps(j,k) = ceil((T1-T0)/tau(j,k));
        tau(j,k) = tau(j,k)-(T0+steps(j,k)*tau(j,k)-T1)/steps(j,k);
        
        % discrete pullback via Interpolation 
        % TODO not implemented
        %pbVm = trace_vertices_RotationMidpoint(NewMesh,-tau(j,k),nsteps_loc, locIntegrate);
        % P_i = assemMat_SemiLag_LFEcheat(NewMesh, pbVm);

        % discrete pullback via L2-projection 
        pbVp = trace_vertices_RotationMidpoint(NewMesh,tau(j,k),nsteps_loc, locIntegrate);
        defMesh = NewMesh;
        defMesh.Coordinates = pbVp(:,[1 2]);
        intersec = aff_elems4(NewMesh,pbVp);
        plot_MeshPatch(NewMesh,defMesh,pbVp,intersec,' ');
        plot_transportMesh(NewMesh,pbVp,'pt ');
        P_p = assemMat_SLpullback_hp(NewMesh, defMesh, Elem2Dof,intersec,quadflag,QuadRule);
               
        % Assemble MASS matrix
        Mex =assemMat_hp(NewMesh,Elem2Dof,@MASS_hp,QuadRule_2D,Shap_2D); 
        M = assemMat_SLmass_hp(NewMesh, defMesh, Elem2Dof,intersec,quadflag,QuadRule);
       
        % inital data
        %u_old = assemCochain_0f(NewMesh,u_Handle,gauleg(0,1,10),0);
        u_old = Mex\assemLoad_hp(NewMesh,Elem2Dof,QuadRule_2D,Shap_2D,u_Handle);
        % Plot hp-FEM solution
        plot_hp(u_old,NewMesh,Elem2Dof,polydegree);
        colorbar; CLim1=get(gca,'CLim');

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
            %u_old = RM\(RM'\(P_p*u_old));
            u_old = M\(P_p*u_old);
            
            %Update vectors
            %H_old=H_new;
            %plot_QFE(u_old,NewMesh);
            %rectangle('Position',[-0.5,-0.5,1,1],'Curvature',[1,1])
            %set(gca,'CLim',CLim1); colorbar;

        end % timestep

        plot_hp(u_old,NewMesh,Elem2Dof,polydegree);
        %rectangle('Position',[-0.5,-0.5,1,1],'Curvature',[1,1])
        set(gca,'CLim',CLim1); colorbar;
       
        error(j,k) = L2Err_hp(NewMesh,u_old,Elem2Dof,QuadRule_2D,Shap_2D,sol_Handle,T1); 
    end %CFLS
end % mesh


data = [h error [zeros(1,size(error,2)); diff(log(error))./diff(log(h*ones(1,size(error,2))))]]

sdata = size(data,1);
fid = fopen(['./results/hp_RotHillscalarL2rate_',locIntegrate,'_',Problem,'.txt'], 'wt');
fprintf(fid,'h');
fprintf(fid, '%4.2d ',CFLs);
fprintf(fid, '\n');
for i = 1:(sdata)
    fprintf(fid, '%12.8e ', data(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n');
fclose(fid)
