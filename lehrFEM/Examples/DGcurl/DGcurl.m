% Run script for discontinuous Galerkin finite element solver
% Numerical results will be saved in the file DGcurl.txt

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  

  % Initialize constants
  fid = fopen('DGcurl.txt','w');
  
  alpha1 = 10;
  alpha2 = 10;
  S = 1;        % Symmetric (+1) or antisymmetric (-1) discretization
  NREFS = 5;     % Number of red refinement steps
  clear Mesh
  SIGMA1 = @(P0,P1,varargin)alpha1/norm(P1-P0);  % interior Edge weight function
  SIGMA2 = @(P0,P1,varargin)alpha2/norm(P1-P0);  % boundary Edge weight function
  fprintf(fid,'DG version with IP-DGFEM with alpha = %10.5f, beta = %10.5f\n',alpha2,alpha1); 
  fprintf(fid,'\n');
  O_Handle = @(x,varargin)0*ones(size(x,1),1);
  for M_flag = 1:2
  for F_flag = [1 3 4]
  clear Mesh 
  if F_flag == 1
      G1 = @(x,varargin)-2*ones(size(x,1),1);  
      G2 = @(x,varargin)2*ones(size(x,1),1);   
      UD1 = @(x,varargin)x(:,2).*(1+x(:,2));   
      UD2 = @(x,varargin)x(:,1).*(1-x(:,1));  
      Curl_U = @(x,varargin)-2*(x(:,1)+x(:,2));
      Div_U = @(x,varargin)0*ones(size(x,1),1);
      fprintf(fid,'Regular function 1 ');
  elseif F_flag == 2
      G1 = @(x,varargin)0*ones(size(x,1),1);  
      G2 = @(x,varargin)0*ones(size(x,1),1);   
      UD1 = @(x,varargin)16/3*(x(:,1).^2+x(:,2).^2).^(13/6).*sin(13*angle(x(:,2),x(:,1))/3);  
      UD2 = @(x,varargin)16/3*(x(:,1).^2+x(:,2).^2).^(13/6).*cos(13*angle(x(:,2),x(:,1))/3); 
      Curl_U = @(x,varargin)0*ones(size(x,1),1);
      Div_U = @(x,varargin)0*ones(size(x,1),1);
      fprintf(fid,'Regular function 2 ');
  elseif F_flag == 3
      G1 = @(x,varargin)0*ones(size(x,1),1);  
      G2 = @(x,varargin)0*ones(size(x,1),1);   
      UD1 = @(x,varargin)4/3*(x(:,1).^2+x(:,2).^2).^(1/6).*sin(angle(x(:,2),x(:,1))/3);  
      UD2 = @(x,varargin)4/3*(x(:,1).^2+x(:,2).^2).^(1/6).*cos(angle(x(:,2),x(:,1))/3);  
      Curl_U = @(x,varargin)0*ones(size(x,1),1);
      Div_U = @(x,varargin)0*ones(size(x,1),1);
      fprintf(fid,'Singular function 1 ');
  elseif F_flag == 4
      G1 = @(x,varargin)0*ones(size(x,1),1);  
      G2 = @(x,varargin)0*ones(size(x,1),1);   
      UD1 = @(x,varargin)2/3*(x(:,1).^2+x(:,2).^2).^(-1/6).*sin(-angle(x(:,2),x(:,1))/3);    
      UD2 = @(x,varargin)2/3*(x(:,1).^2+x(:,2).^2).^(-1/6).*cos(-angle(x(:,2),x(:,1))/3);  
      Curl_U = @(x,varargin)0*ones(size(x,1),1);
      Div_U = @(x,varargin)0*ones(size(x,1),1);
      fprintf(fid,'Singular function 2 ');
  elseif F_flag == 5 
      G1 = @(x,varargin)0*ones(size(x,1),1);
      G2 = @(x,varargin)0*ones(size(x,1),1);
      UD1 = @(x,varargin)-x(:,1).^2+x(:,2).^2-2*x(:,1).*x(:,2);
      UD2 = @(x,varargin)-x(:,1).^2+x(:,2).^2+2*x(:,1).*x(:,2);
      fprintf(fid,'Regular function 4 ');
  else    
      G1 = @(x,varargin)-2*exp(x(:,1)).*sin(x(:,2)); 
      G2 = @(x,varargin)-2*exp(x(:,1)).*cos(x(:,2));   
      UD1 = @(x,varargin)-exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2)));  
      UD2 = @(x,varargin)exp(x(:,1)).*x(:,2).*sin(x(:,2));
      Curl_U = @(x,varargin)2*exp(x(:,1)).*cos(x(:,2)); 
      Div_U = @(x,varargin)0*ones(size(x,1),1);
      fprintf(fid,'Regular function 3 ');
  end
  
  % Initialize mesh
  if M_flag == 1 
  Mesh.Coordinates = [0 0;1 0; 1 1; 0 1];
  Mesh.Elements = [1 2 3; 1 3 4];
  fprintf(fid,'in square domain\n');
  elseif M_flag == 2
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  fprintf(fid,'in L-shaped domain\n');
  else
  load LEBLmesh.mat 
  Mesh = Mesh5;
  Mesh = rmfield(Mesh,'BdFlags');
  fprintf(fid,'in Graded L-shaped domain\n');   
  end
  
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh); 
  
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  if M_flag~=3 & M_flag~=4
  Mesh = refine_REG(Mesh);
  Mesh = refine_REG(Mesh);
  end
  
  fprintf(fid,'Mesh Width   L2 Norm of solution   Relative L2 Error     Order     DG Norm of solution   Relative DG norm of Error     Order\n');
  fprintf(fid,'==========   ===================   =================   =========   ===================   =========================   =========\n');
  
  for i = 1:NREFS
  
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  
  
  % Assemble matrices and load vectors (discontinuous Lagrangian elements)
  
  QuadRule_1D = gauleg(0,1,10);
  QuadRule_2D = P3O3();
  
  [IC1,JC1,Avol_Curl] = assemMat_Vol_DG2(Mesh,@STIMA_Curl_LFE2);
  [ID1,JD1,Avol_Div] = assemMat_Vol_DG2(Mesh,@STIMA_Div_LFE2);
  
  [It2,Jt2,Jinn_t] = assemMat_Inn_DG2(Mesh,@STIMA_InnPen_Tangent_DGLFE2,SIGMA2);
  [In2,Jn2,Jinn_n] = assemMat_Inn_DG2(Mesh,@STIMA_InnPen_Normal_DGLFE2,SIGMA1);
  [IC2,JC2,Ainn_Curl] = assemMat_Inn_DG2(Mesh,@STIMA_Curl_Inn_DGLFE2,S);   
  [ID2,JD2,Ainn_Div] = assemMat_Inn_DG2(Mesh,@STIMA_Div_Inn_DGLFE2,S);
    
  [It3,Jt3,Jbnd] = assemMat_Bnd_DG2(Mesh,@STIMA_BndPen_Tangent_DGLFE2,SIGMA2); 
  [IC3,JC3,Abnd_Curl] = assemMat_Bnd_DG2(Mesh,@STIMA_Curl_Bnd_DGLFE2,S);   
 
  Lvol = assemLoad_Vol_DG2(Mesh,@LOAD_Vol_DGLFE2,QuadRule_2D,G1,G2);
  Lbnd = assemLoad_Bnd_DG2(Mesh,@LOAD_Curl_Bnd_Tangent_DGLFE2,QuadRule_1D,S,SIGMA2,UD1,UD2);
  h=get_MeshWidth(Mesh); 
  
  % Create system matrix
       
  A = sparse([JC1; JD1; JC2; JD2; JC3; Jt2; Jn2; Jt3], ...
             [IC1; ID1; IC2; ID2; IC3; It2; In2; It3],...
             [Avol_Curl; Avol_Div; Ainn_Curl; Ainn_Div; Abnd_Curl; Jinn_t; Jinn_n; Jbnd]);
  L = Lvol + Lbnd;
 
  % Solve the linear system
  U = A\L;
  
  M_Curl = sparse([JC1], ...
             [IC1], ...
             [Avol_Curl]);    
  
  M_Div = sparse([JD1], ...
             [ID1], ...
             [Avol_Div]);    
  
  X = sparse([Jt2;Jn2;Jt3], ...
             [It2;Jn2;It3], ...
             [1/alpha2*Jinn_t;1/alpha1*Jinn_n;0/alpha2*Jbnd]);    
  Y = sparse([Jn2;Jt3], ...
             [In2;It3], ...
             [0/alpha2*Jinn_t;1/alpha2*Jbnd]);   
             
  % Interpolation of exaction solution           
  UE = L2_interp_DGLFE(Mesh,P7O6(),UD1,UD2);
  
  % Interpolation of exaction solution
  
  if i == NREFS+1
  plot_DGLFE(U(1:2:end),Mesh);
  colorbar;

  plot_DGLFE(U(2:2:end),Mesh);  
  colorbar;
  end
  
  % Compute the errors
  
  M_W(i) = get_MeshWidth(Mesh);
  Norm(i) =  L2Err_DGLFE2(Mesh,U,P7O6(),O_Handle,O_Handle);
  L2Err(i) = L2Err_DGLFE2(Mesh,U,P7O6(),UD1,UD2);
  DGNorm(i) = EnergyErr_DGcurl(Mesh,U,P7O6,O_Handle,O_Handle,Curl_U,Div_U,SIGMA1,SIGMA2,alpha1,alpha2);
  DGErr(i) = EnergyErr_DGcurl(Mesh,U,P7O6,UD1,UD2,Curl_U,Div_U,SIGMA1,SIGMA2,alpha1,alpha2);
  
  if i==1
      fprintf(fid,'%8.4f & %14.4e      &%17.4e    &     -     & %17.4e   &    %17.4e  &        -   \\\\ \n',M_W(i),Norm(i),L2Err(i),DGNorm(i),DGErr(i));
  else
      fprintf(fid,'%8.4f & %14.4e      &%17.4e    &  %8.5f & %17.4e   &    %17.4e  &     %8.5f\\\\ \n',M_W(i),Norm(i),L2Err(i),log(L2Err(i-1)/L2Err(i))/log(2),DGNorm(i),DGErr(i),log(DGErr(i-1)/DGErr(i))/log(2));
  end
  
  save Soln_Curl.mat U Mesh
  Mesh = refine_REG(Mesh);
  clear A L U
  end
  
    if M_flag==1 
      if F_flag==1
          save DGcurl11.mat M_W L2Err
      elseif F_flag==3
          save DGcurl13.mat M_W L2Err
      elseif F_flag==4
          save DGcurl14.mat M_W L2Err
      end
    else
      if F_flag==1
          save DGcurl21.mat M_W L2Err
      elseif F_flag==3
          save DGcurl23.mat M_W L2Err
      elseif F_flag==4
          save DGcurl24.mat M_W L2Err
      end
    end
    
  fprintf(fid,'\n');
  end
  end
  
  fclose(fid);

  % Clear memory
    
  clear all
