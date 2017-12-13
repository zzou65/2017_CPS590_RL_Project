% test bilinear forms
clear Mesh

NREFS=7;
JIG=2;
d=getData(13);
QuadRule=Duffy(TProd(gauleg(0,1,10)));

% Initialize mesh
MU_HANDLE=@(x,varargin)1;
W0_handle=@(x,varargin)x(:,2)+x(:,1);
%W_handle=@(x,varargin)[(1+x(:,1)).^2 5*x(:,2).*x(:,1)];
W_handle=@(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];

result=innerproduct(d.SOL2_Handle,W_handle,0,1,0,1)

Mesh.Coordinates =[0 0; 1 0; 1 1; 0 1];
Mesh.Elements = [1 2 4;2 3 4];
   % Mesh = load_Mesh('Coordinates.dat','Elements.dat');
    
Mesh = add_Edges(Mesh);
Mesh = add_Edge2Elem(Mesh);
    
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = [-1 -1 -2 -2];
    
    for i = 1:NREFS
    
        Mesh=refine_REG(Mesh);
        Mesh=add_Edge2Elem(Mesh);
        Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
        
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
          Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
        end
        
        M_W1F = assemMat_W1F(New_Mesh,@MASS_W1F,MU_HANDLE, P3O3());
        %MassOne=assemMat_Mass1fD(Mesh);                         % diagonal Mass Matrix of 1-Forms 
        M_LFE=assemMat_LFE(Mesh,@MASS_LFE);
        GradContr=assemMat_W1F(New_Mesh,@STIMA_GradContr, d.V_Handle, gauleg(0,1,10));
        %GradContr2=assemMat_W1F(New_Mesh,@STIMA_GradContr2, d.V_Handle, QuadRule);
        %SUPG=assemMat_W1F(New_Mesh,@STIMA_SUPG_W1F, d.V_Handle, QuadRule);
        %Lie=assemMat_W1F(New_Mesh,@STIMA_Lie_W1F, d.V_Handle);
        %ContrRot=assemMat_W1F(New_Mesh,@STIMA_ContrRot, d.V_Handle, QuadRule);
        %ContrOne=assemMat_Contr1f(New_Mesh,d.V_Handle);
        TopGrad=assemMat_TopGrad(New_Mesh);   % topological Gradient
        ContrOne=assemMat_Contr1f(New_Mesh,d.V_Handle);  % contraction of one forms 
        
        cVUBd=assemBndCochain_0f(New_Mesh,[-2 -3],d.V_U_EX_Handle); % boundary data on outflow boundary
       
        L1=assemLoad_W1F(New_Mesh,P7O6(),d.U_EX_Handle);
            
U1=M_W1F\L1;
        L2=assemLoad_W1F(New_Mesh,P7O6(),W_handle);
        U2=M_W1F\L2;
        [U2'*GradContr*U1 U2'*(M_W1F*TopGrad*cVUBd+M_W1F*TopGrad*ContrOne*U1) result]
        %err(i)=abs(U2'*GradContr*U1-result);
        err(i)=abs(U2'*(M_W1F*TopGrad*cVUBd+M_W1F*TopGrad*ContrOne*U1)-result);
        %err(i)=L2Err_W1F_mod(New_Mesh,(M_W1F*TopGrad*cVUBd+M_W1F*TopGrad*ContrOne*U1),QuadRule,d.SOL2_Handle);
        h(i)=get_MeshWidth(New_Mesh);
        plot_Norm_W1F(M_W1F*TopGrad*cVUBd+M_W1F*TopGrad*ContrOne*U1, Mesh); colorbar;
        % figure; plot_W1F(GradContr*U1, Mesh)
        
    end
    fig = figure('Name','rate');
    plot(h,err); grid('on');
    set(gca,'XScale','log','YScale','log');
    xlabel('{\bf h}');
    ylabel('{\bf Error}');
    p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS)),1);
    add_Slope(gca,'SouthWest',p(1));
    