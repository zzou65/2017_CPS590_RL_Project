% Run script for discontinuous Galerkin finite element solver
% Numerical results will be saved in the file DGlapl.txt

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  

  % Initialize constants
 
  fid = fopen('DGlapl.txt','w');

  alpha = 10;
  S = 1;        % Symmetric (+1) or antisymmetric (-1) discretization
  NREFS = 5;     % Number of red refinement steps
  fprintf(fid,'DG Laplacian version with IP-DGFEM\n');  
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
      grad11 = @(x,varargin)0*ones(size(x,1),1);
      grad12 = @(x,varargin)1+2*x(:,2);
      grad21 = @(x,varargin)1-2*x(:,1);
      grad22 = grad11;
      fprintf(fid,'Regular function 1 ');
  elseif F_flag == 2
      G1 = @(x,varargin)0*ones(size(x,1),1);  
      G2 = @(x,varargin)0*ones(size(x,1),1);   
      UD1 = @(x,varargin)16/3*(x(:,1).^2+x(:,2).^2).^(13/6).*sin(13*angle(x(:,2),x(:,1))/3);  
      UD2 = @(x,varargin)16/3*(x(:,1).^2+x(:,2).^2).^(13/6).*cos(13*angle(x(:,2),x(:,1))/3);
      grad11 = @(x,varargin)16*13/9*(x(:,1).^2+x(:,2).^2).^(10/6).*sin(10*angle(x(:,2),x(:,1))/3); 
      grad12 = @(x,varargin)16*13/9*(x(:,1).^2+x(:,2).^2).^(10/6).*cos(10*angle(x(:,2),x(:,1))/3); 
      grad21 = @(x,varargin)16*13/9*(x(:,1).^2+x(:,2).^2).^(10/6).*cos(10*angle(x(:,2),x(:,1))/3); 
      grad22 = @(x,varargin)-16*13/9*(x(:,1).^2+x(:,2).^2).^(10/6).*sin(10*angle(x(:,2),x(:,1))/3); 
      fprintf(fid,'Regular function 2  ');
  elseif F_flag == 3
      G1 = @(x,varargin)0*ones(size(x,1),1);  
      G2 = @(x,varargin)0*ones(size(x,1),1);   
      UD1 = @(x,varargin)4/3*(x(:,1).^2+x(:,2).^2).^(1/6).*sin(angle(x(:,2),x(:,1))/3);  
      UD2 = @(x,varargin)4/3*(x(:,1).^2+x(:,2).^2).^(1/6).*cos(angle(x(:,2),x(:,1))/3);  
      grad11 = @(x,varargin)-4/9*(x(:,1).^2+x(:,2).^2).^(-2/6).*sin(2*angle(x(:,2),x(:,1))/3); 
      grad12 = @(x,varargin)4/9*(x(:,1).^2+x(:,2).^2).^(-2/6).*cos(2*angle(x(:,2),x(:,1))/3); 
      grad21 = @(x,varargin)4/9*(x(:,1).^2+x(:,2).^2).^(-2/6).*cos(2*angle(x(:,2),x(:,1))/3); 
      grad22 = @(x,varargin)4/9*(x(:,1).^2+x(:,2).^2).^(-2/6).*sin(2*angle(x(:,2),x(:,1))/3); 
      fprintf(fid,'Singular function 1 ');
  elseif F_flag == 4
      G1 = @(x,varargin)0*ones(size(x,1),1);  
      G2 = @(x,varargin)0*ones(size(x,1),1);   
      UD1 = @(x,varargin)2/3*(x(:,1).^2+x(:,2).^2).^(-1/6).*sin(-angle(x(:,2),x(:,1))/3);    
      UD2 = @(x,varargin)2/3*(x(:,1).^2+x(:,2).^2).^(-1/6).*cos(-angle(x(:,2),x(:,1))/3); 
      grad11 = @(x,varargin)2/9*(x(:,1).^2+x(:,2).^2).^(-4/6).*sin(4*angle(x(:,2),x(:,1))/3); 
      grad12 = @(x,varargin)-2/9*(x(:,1).^2+x(:,2).^2).^(-4/6).*cos(4*angle(x(:,2),x(:,1))/3); 
      grad21 = @(x,varargin)-2/9*(x(:,1).^2+x(:,2).^2).^(-4/6).*cos(4*angle(x(:,2),x(:,1))/3); 
      grad22 = @(x,varargin)-2/9*(x(:,1).^2+x(:,2).^2).^(-4/6).*sin(4*angle(x(:,2),x(:,1))/3); 
      %UD1 = @(x,varargin)sing_fcn2(x);
      %UD2 = @(x,varargin)sing_fcn(x);
      fprintf(fid,'Singular function 2 ');
  elseif F_flag == 5 
      G1 = @(x,varargin)0*ones(size(x,1),1);
      G2 = @(x,varargin)0*ones(size(x,1),1);
      UD1 = @(x,varargin)-x(:,1).^2+x(:,2).^2-2*x(:,1).*x(:,2);
      UD2 = @(x,varargin)-x(:,1).^2+x(:,2).^2+2*x(:,1).*x(:,2);
      grad11 = @(x,varargin)-2*x(:,1)-2*x(:,2);
      grad12 = @(x,varargin)-2*x(:,1)+2*x(:,2);
      grad21 = @(x,varargin)-2*x(:,1)+2*x(:,2);
      grad22 = @(x,varargin)2*x(:,1)+2*x(:,2);
      fprintf(fid,'Regular function 4');
  elseif F_flag ==6
      G1 = @(x,varargin)-2*exp(x(:,1)).*sin(x(:,2)); 
      G2 = @(x,varargin)-2*exp(x(:,1)).*cos(x(:,2));   
      UD1 = @(x,varargin)-exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2)));  
      UD2 = @(x,varargin)exp(x(:,1)).*x(:,2).*sin(x(:,2));
      grad11 = @(x,varargin)-exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2)));
      grad12 = @(x,varargin)-exp(x(:,1)).*(2*cos(x(:,2))-x(:,2).*sin(x(:,2)));
      grad21 = @(x,varargin)exp(x(:,1)).*x(:,2).*sin(x(:,2));
      grad22 = @(x,varargin)exp(x(:,1)).*(x(:,2).*cos(x(:,2))+sin(x(:,2)));
      fprintf(fid,'Regular function 3 ');
  else
      G1 = @(x,varargin)-3*sin(x(:,1)).*exp(2*x(:,2));
      G2 = @(x,varargin)-3*cos(x(:,2)).*exp(2*x(:,1));
      UD1 = @(x,varargin)sin(x(:,1)).*exp(2*x(:,2));
      UD2 = @(x,varargin)cos(x(:,2)).*exp(2*x(:,1));
      grad11 = @(x,varargin)cos(x(:,1)).*exp(2*x(:,2));
      grad12 = @(x,varargin)2*sin(x(:,1)).*exp(2*x(:,2));
      grad21 = @(x,varargin)2*cos(x(:,2)).*exp(2*x(:,1));
      grad22 = @(x,varargin)-sin(x(:,2)).*exp(2*x(:,1));
      div_u = @(x,varargin)cos(x(:,1)).*exp(2*x(:,2))-sin(x(:,2)).*exp(2*x(:,1));
      fprintf(fid,'Regular function 3 ');
  end
  
  
                                 
  SIGMA = @(P0,P1,varargin)alpha/norm(P1-P0);  % Edge weight function
  
  % Initialize mesh
  if M_flag == 1 
  Mesh.Coordinates = [0 0;1 0; 1 1; 0 1];
  Mesh.Elements = [1 2 3; 1 3 4];

  fprintf(fid,'in square domain\n');
  elseif M_flag == 2
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  fprintf(fid,'in L-shaped domain\n');
  elseif M_flag == 3
  load gradedL_43.mat Mesh
  fprintf(fid,'in Graded L-shaped domain\n');   
  else
      load tri.mat
      fprintf(fid,'in triangular domain\n');
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
  
  fprintf(fid,'Mesh Width   L2 Norm of exact solution   Relative L2 Error     Order     DG Norm of solution   Relative DG norm of Error     Order\n');
  fprintf(fid,'==========   =========================   =================   =========   ===================   =========================   =========\n');
  
  for i = 1:NREFS
  
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);

  % Assemble matrices and load vectors (discontinuous Lagrangian elements)
  
  QuadRule_1D = gauleg(0,1,10);
  QuadRule_2D = P3O3();
  
  [I1,J1,Avol] = assemMat_Vol_DG2(Mesh,@STIMA_Lapl_Vol_DGLFE2);
  
  [I2,J2,Jinn] = assemMat_Inn_DG2(Mesh,@STIMA_InnPen_DGLFE2,SIGMA);
  [I2,J2,Ainn] = assemMat_Inn_DG2(Mesh,@STIMA_Grad_Inn_DGLFE2,S);
    
  [I3,J3,Jbnd] = assemMat_Bnd_DG2(Mesh,@STIMA_BndPen_Tangent_DGLFE2,SIGMA);  
  [I3,J3,Abnd] = assemMat_Bnd_DG2(Mesh,@STIMA_Grad_Bnd_Tangent_DGLFE2,S);
 
  Lvol = assemLoad_Vol_DG2(Mesh,@LOAD_Vol_DGLFE2,QuadRule_2D,G1,G2);
  Lbnd = assemLoad_Bnd_DG2(Mesh,@LOAD_Grad_Bnd_Tangent_DGLFE2,QuadRule_1D,S,SIGMA,UD1,UD2);
  Lbnd_Du = assemLoad_Bnd_DG2(Mesh,@LOAD_Du_Bnd_Normal_DGLFE2,QuadRule_1D,S,SIGMA,grad11,grad12,grad21,grad22);
  
  % Create system matrix 
 
  A = sparse([J1; J2; J3; J2; J3], ...
             [I1; I2; I3; I2; I3], ...
             [Avol; Ainn; Abnd; Jinn; Jbnd]);
  L = Lvol + Lbnd - Lbnd_Du;
 
  % Solve the linear system
  
  U = A\L;
  
  if i == NREFS+1
  plot_DGLFE(U(1:2:end),Mesh);
  colorbar;

  plot_DGLFE(U(2:2:end),Mesh);  
  colorbar;
  end
  
  % Interpolation of exaction solution
  
  UE = L2_interp_DGLFE(Mesh,P7O6(),UD1,UD2);
  
  % Compute the errors
  L2Err1 = L2Err_DGLFE(Mesh,U(1:2:end),P7O6(),UD1);
  L2Err2 = L2Err_DGLFE(Mesh,U(2:2:end),P7O6(),UD2);
  norm_Ex1 = L2Err_DGLFE(Mesh,zeros(size(U,1)/2,1),P7O6(),UD1);
  norm_FE1 = L2Err_DGLFE(Mesh,U(1:2:end),P7O6(),O_Handle);
  norm_Ex2 = L2Err_DGLFE(Mesh,zeros(size(U,1)/2,1),P7O6(),UD2);
  norm_FE2 = L2Err_DGLFE(Mesh,U(2:2:end),P7O6(),O_Handle);
  M_W(i) = get_MeshWidth(Mesh);
  Norm_Ex(i) = sqrt(norm_Ex1^2+norm_Ex2^2);
  Norm(i) = sqrt(norm_FE1^2+norm_FE2^2);
  L2Err(i) = sqrt(L2Err1^2+L2Err2^2);
  A1 = sparse([J1;J3], ...
             [I1;I3], ...
             [Avol;0*Jbnd]);
  A2 = sparse([J2], ...
             [I2], ...
             [Jinn]); 
  A3 = sparse([J3;J1], ...
             [I3;I1], ...
             [Jbnd;0*Avol]);        
  DGNorm(i) = sqrt(U'*A1*U+U'*A2*U+U'*A3*U);
  DGErr(i) = sqrt((U-UE)'*A1*(U-UE)+(U-UE)'*A2*(U-UE)+(U-UE)'*A3*(U-UE));
  
  if i==1
      fprintf(fid,'%8.4f & %20.4e   &   %17.4e    &     -     & %17.4e   &    %17.4e  &        -   \\\\ \n',M_W(i),Norm(i),L2Err(i),DGNorm(i),DGErr(i));
  else
      fprintf(fid,'%8.4f & %20.4e   &   %17.4e    &  %8.5f & %17.4e   &    %17.4e  &     %8.5f\\\\ \n',M_W(i),Norm(i),L2Err(i),log(L2Err(i-1)/L2Err(i))/log(2),DGNorm(i),DGErr(i),log(DGErr(i-1)/DGErr(i))/log(2));
  end
  Mesh = refine_REG(Mesh);
  clear A L U
  end
  fprintf(fid,'\n');
  end
  end

  fclose(fid);
  
  % Clear memory
    
  clear all
