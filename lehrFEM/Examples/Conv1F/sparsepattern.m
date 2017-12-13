% Run script discrete differential forms.

%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    %close all;
 %   clear all;
    
    % Initialize constant
    %e_a* (v x curl u)+e_b*(grad v.u)+e_c* curlcurl u +e_m*u
    e_a=1;                            %(v x curl u)   
    e_b=0;                            %(grad v.u)
    e_c=10^-4;                      %curlcurl u 
    e_m=1;                            %u
    %d=getData(12);                   % struct containing Coordinates, Elements and various function handles
    EPSI_Handle = @(x,varargin)ones(size(x,1),1);
    MU_HANDLE = @(x,varargin)1;
    QuadRule = Duffy(TProd(gauleg(0,1,10)));
   
    JIG =2;  
    
    % Initialize mesh
%   
    NREFS =3;
    Mesh.Coordinates = [-1.5 -1.5; 1.5 -1.5; 1.5 1.5; -1.5 1.5];
    Mesh.Elements = [1 2 3;1 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = [-1 -1 -2 -2];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    for i = 1:NREFS
    i
        Mesh=refine_REG(Mesh);
        Mesh=add_Edge2Elem(Mesh);
        
        % Mesh preprocessing
                
   switch(JIG)
        case 1
          New_Mesh = Mesh;      
        case 2
          Loc = get_BdEdges(Mesh);
          Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
          FixedPos = zeros(size(Mesh.Coordinates,1),1);
          FixedPos(Loc) = 1;
          New_Mesh = jiggle(Mesh,FixedPos);   
        case 3
          Loc = get_BdEdges(Mesh);
          Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
          FixedPos = zeros(size(Mesh.Coordinates,1),1);
          FixedPos(Loc) = 1;
          New_Mesh = smooth(Mesh,FixedPos);
   end

   TopGrad=assemMat_TopGrad(New_Mesh);   % topological Gradient
   TopRot=assemMat_TopRot(New_Mesh);      % topological Rotation
   
   MassZero=assemMat_Mass0fD(Mesh);         %diagonal 
   
    % several mass matrices
 
    M_W1F = assemMat_W1F(New_Mesh,@MASS_W1F,MU_HANDLE, P3O3());  % Mass Matrix Whitney-1-Forms
    MassOne=assemMat_Mass1fD(Mesh);  
  %  REG=M_W1F*TopGrad*MassZero*TopGrad'*M_W1F;
    REG=MassOne*TopGrad*MassZero*TopGrad'*MassOne;
    
    end
    h=Thetas;
    plot_Mesh(New_Mesh)
    fig = figure('Name','Discretization error');
    plot(h,abs(err1),'ro-',h,abs(err2),'go-',h,abs(err3),'bo-'); grid('on');
    set(gca,'YScale','log');
    xlabel('{\bf 2\pi/ \Theta}');
    ylabel('{\bf norms of inverse}');
    
    legend('standard','supg','upwind','Location','NorthEast') 
clear all