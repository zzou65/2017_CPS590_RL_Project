% Run script for discontinuous Galerkin finite element solver
% Numerical results will be saved in the file MIXDG.txt

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  

  % Initialize constants
  fid = fopen('MIXDG.txt','w');
  
  alpha1 = 1;   % parameter beta
  alpha2 = 10;  % parameter alpha
  S = 1;        % Symmetric (+1) or antisymmetric (-1) discretization
  NREFS = 5;     % Number of red refinement steps
  clear Mesh
  fprintf(fid,'Mixed DG with alpha = %10.5f, beta = %10.5f\n',alpha2,alpha1);  
  fprintf(fid,'\n');
  SIGMA1 = @(P0,P1,varargin)alpha1*norm(P1-P0);  % interior Edge weight function
  SIGMA2 = @(P0,P1,varargin)alpha2/norm(P1-P0);  % boundary Edge weight function
  
  O_Handle = @(x,varargin)0*ones(size(x,1),1);   % zero function
  O_Handle2 = @(x,varargin)0*[ones(size(x,1),1) ones(size(x,1),1)];
  
  for M_flag = 1:2
  for F_flag = [1 3 4]
  clear Mesh 
  
  % Initialize source function and boundary datum for different exact solution
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
      G1 = @(x,varargin)-2*exp(x(:,1)).*sin(x(:,2))+pi*2*cos(pi*(x(:,1)-1)*2).*sin(pi*(x(:,2)-1)*2);  
      G2 = @(x,varargin)-2*exp(x(:,1)).*cos(x(:,2))+pi*2*sin(pi*(x(:,1)-1)*2).*cos(pi*(x(:,2)-1)*2);    
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
  
  % Start computation 
  for i = 1:NREFS 
  
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  
  % Assemble matrices and load vectors (discontinuous Lagrangian elements)
  QuadRule_1D = gauleg(0,1,10);
  QuadRule_2D = P3O3();
  
  [I1,J1,Avol_Curl] = assemMat_Vol_DG2(Mesh,@STIMA_Curl_LFE2);
  [I1,J1,Avol_Div] = assemMat_Vol_DG2(Mesh,@STIMA_Div_LFE2);
  
  [IB,JB,B] = assemMat_MIXDG(Mesh,@STIMA_Grad_MIXDG);
  D = assemMat_LFE(Mesh,@MASS_Lump_Inv_LFE);
  B = sparse([IB],[JB],[B]);
  
  nCoordinates = size(Mesh.Coordinates,1);
  Loc = get_BdEdges(Mesh);
  DEdges = Loc(Mesh.BdFlags(Loc) == -1);
  DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
  FreeDofs = setdiff(1:nCoordinates,DNodes);
  
  B = B(:,FreeDofs);
  D = D(FreeDofs,FreeDofs);
  
  T = B*D*transpose(B);
  [IT,JT,T] = find(T);
  
  [I2,J2,Jinn_n] = assemMat_Inn_DG2(Mesh,@STIMA_InnPen_Normal_DGLFE2,SIGMA1);
  [I2,J2,Jinn_t] = assemMat_Inn_DG2(Mesh,@STIMA_InnPen_Tangent_DGLFE2,SIGMA2);
  
  [I2,J2,Ainn_Curl] = assemMat_Inn_DG2(Mesh,@STIMA_Curl_Inn_DGLFE2,S);   
  
  [I3,J3,Jbnd] = assemMat_Bnd_DG2(Mesh,@STIMA_BndPen_Tangent_DGLFE2,SIGMA2); 
  [I3,J3,Abnd_Curl] = assemMat_Bnd_DG2(Mesh,@STIMA_Curl_Bnd_DGLFE2,S);   
 
  Lvol = assemLoad_Vol_DG2(Mesh,@LOAD_Vol_DGLFE2,QuadRule_2D,G1,G2);
  Lbnd = assemLoad_Bnd_DG2(Mesh,@LOAD_Curl_Bnd_Tangent_DGLFE2,QuadRule_1D,S,SIGMA2,UD1,UD2);
  
  h=get_MeshWidth(Mesh);
  
  % Create system matrix
       
  A = sparse([J1; J2; J3; JT], ...
             [I1; I2; I3; IT], ...
             [Avol_Curl+alpha1*h^2*Avol_Div; Ainn_Curl+Jinn_t+Jinn_n; Abnd_Curl+Jbnd; T]);
  L = Lvol+Lbnd;
  
  U = A\L;
  
  % Plot the numerical solution in the last refinement
  
  if i == NREFS+1
    plot_DGLFE(U(1:2:end),Mesh);
    colorbar;
    plot_DGLFE(U(2:2:end),Mesh);  
    colorbar;
  end

  % Compute and print the errors

  M_W(i) = get_MeshWidth(Mesh);
  Norm(i) = L2Err_DGLFE2(Mesh,U,P7O6(),O_Handle,O_Handle);
  L2Err(i) = L2Err_DGLFE2(Mesh,U,P7O6(),UD1,UD2);
  DGNorm(i) = EnergyErr_MIXDG(Mesh,U,P7O6,O_Handle,O_Handle,Curl_U,Div_U,SIGMA1,SIGMA2,alpha1,alpha2);
  DGErr(i) = EnergyErr_MIXDG(Mesh,U,P7O6,UD1,UD2,Curl_U,Div_U,SIGMA1,SIGMA2,alpha1,alpha2);

  % If one want to print the errors of muliplier p, use the following
  %
  %nv = size(Mesh.Coordinates,1);
  %p = zeros(nv,1);
  %p(FreeDofs) = D*B'*U;
  %L2Err(i) = L2Err_LFE(Mesh,p,P7O6(),O_Handle);
  %Norm(i) = L2Err(i);
  %DGErr(i) = H1Err_LFE(Mesh,p,P7O6(),O_Handle,O_Handle2);
  %DGNorm(i) = DGErr(i);

  if i==1
      fprintf(fid,'%8.4f & %14.4e   &   %17.4e    &     -     & %17.4e   &    %17.4e  &        -   \\\\ \n',M_W(i),Norm(i),L2Err(i),DGNorm(i),DGErr(i));
  else
      fprintf(fid,'%8.4f & %14.4e   &   %17.4e    &  %8.5f & %17.4e   &    %17.4e  &     %8.5f\\\\ \n',M_W(i),Norm(i),L2Err(i),log(L2Err(i-1)/L2Err(i))/log(2),DGNorm(i),DGErr(i),log(DGErr(i-1)/DGErr(i))/log(2));
  end
  
  % Refine the mesh after each computation
  Mesh = refine_REG(Mesh);
  clear A L U
  end
  
    % Save the errors and mesh width for plotting purpose 
    if M_flag==1 
      if F_flag==1
          save MIXDG11.mat M_W L2Err
      elseif F_flag==3
          save MIXDG13.mat M_W L2Err
      elseif F_flag==4
          save MIXDG14.mat M_W L2Err
      end
    else
      if F_flag==1
          save MIXDG21.mat M_W L2Err
      elseif F_flag==3
          save MIXDG23.mat M_W L2Err
      elseif F_flag==4
          save MIXDG24.mat M_W L2Err
      end
    end
    
  fprintf(fid,'\n');
  end
  end
  
  fclose(fid);

  % Clear memory
    
  clear all
