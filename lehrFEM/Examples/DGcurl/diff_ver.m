% Run script for comparing the L2 difference of solution of the two methods
% in appendix
% Numerical results will be saved in the file diff_ver.txt

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  fid = fopen('diff_ver.txt','w');
  
  s_div = 1;
  alpha1 = 10;
  alpha2 = 10;
  S = -1;        % Symmetric (+1) or antisymmetric (-1) discretization
  NREFS = 5;     % Number of red refinement steps
  clear Mesh
  O_Handle = @(x,varargin)0*ones(size(x,1),1);
  fprintf(fid,'Difference of solutions of Laplacian and Curl Curl version with IP-DGFEM\n');  
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
      grad12 = grad11;
      grad21 = grad11;
      grad22 = grad11;
      fprintf(fid,'Regular function 1 ');
  elseif F_flag == 2
      G1 = @(x,varargin)-2*ones(size(x,1),1);  
      G2 = @(x,varargin)2*ones(size(x,1),1);   
      UD1 = @(x,varargin)x(:,2).*(2+x(:,2));   
      UD2 = @(x,varargin)x(:,1).*(2-x(:,1));  
      fprintf(fid,'Regular function 2 ');
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
      fprintf(fid,'Singular function 2 ');
 elseif F_flag == 5 
      G1 = @(x,varargin)D1(norm(x)).*sin(angle(x(:,2),x(:,1))/3+pi).*[norm(x)>.25 & norm(x)<.75]-D2(norm(x)).*sin(7*angle(x(:,2),x(:,1))/3+pi).*[norm(x)>.25 & norm(x)<.75];
      G2 = @(x,varargin)D1(norm(x)).*cos(angle(x(:,2),x(:,1))/3+pi).*[norm(x)>.25 & norm(x)<.75]+D2(norm(x)).*cos(7*angle(x(:,2),x(:,1))/3+pi).*[norm(x)>.25 & norm(x)<.75];
      UD1 = @(x,varargin)-4/3*(x(:,1).^2+x(:,2).^2).^(1/6).*sin(angle(x(:,2),x(:,1))/3+pi).*[norm(x)<=.25]-4/3*(x(:,1).^2+x(:,2).^2).^(1/6).*sin(angle(x(:,2),x(:,1))/3+pi).*phi(norm(x)).*[norm(x)>.25 & norm(x)<.75]+norm(x).^(4/3).*phi1(norm(x)).*cos(4*angle(x(:,2),x(:,1))/3+pi).*sin(angle(x(:,2),x(:,1))).*[norm(x)>.25 & norm(x)<.75];    
      UD2 = @(x,varargin)-4/3*(x(:,1).^2+x(:,2).^2).^(1/6).*cos(angle(x(:,2),x(:,1))/3+pi).*[norm(x)<=.25]-4/3*(x(:,1).^2+x(:,2).^2).^(1/6).*cos(angle(x(:,2),x(:,1))/3+pi).*phi(norm(x)).*[norm(x)>.25 & norm(x)<.75]-norm(x).^(4/3).*phi1(norm(x)).*cos(4*angle(x(:,2),x(:,1))/3+pi).*cos(angle(x(:,2),x(:,1))).*[norm(x)>.25 & norm(x)<.75];    
      fprintf(fid,'Singular function 3 ');
  else    
      G1 = @(x,varargin)-G11(norm(x)).*sin(angle(x(:,2),x(:,1))/3+pi).*[norm(x)>.25 & norm(x)<.75]-G12(norm(x)).*sin(5*angle(x(:,2),x(:,1))/3+pi).*[norm(x)>.25 & norm(x)<.75];
      G2 = @(x,varargin)-G21(norm(x)).*cos(angle(x(:,2),x(:,1))/3+pi).*[norm(x)>.25 & norm(x)<.75]-G22(norm(x)).*cos(5*angle(x(:,2),x(:,1))/3+pi).*[norm(x)>.25 & norm(x)<.75];
      UD1 = @(x,varargin)2/3*(x(:,1).^2+x(:,2).^2).^(-1/6).*sin(angle(x(:,2),x(:,1))/3+pi).*[norm(x)<=.25]+2/3*(x(:,1).^2+x(:,2).^2).^(-1/6).*sin(angle(x(:,2),x(:,1))/+pi).*phi(norm(x)).*[norm(x)>.25 & norm(x)<.75]+norm(x).^(2/3).*phi1(norm(x)).*cos(2*angle(x(:,2),x(:,1))/3+pi).*sin(angle(x(:,2),x(:,1))).*[norm(x)>.25 & norm(x)<.75];    
      UD2 = @(x,varargin)-2/3*(x(:,1).^2+x(:,2).^2).^(-1/6).*cos(angle(x(:,2),x(:,1))/3+pi).*[norm(x)<=.25]-2/3*(x(:,1).^2+x(:,2).^2).^(-1/6).*cos(angle(x(:,2),x(:,1))/3+pi).*phi(norm(x)).*[norm(x)>.25 & norm(x)<.75]-norm(x).^(2/3).*phi1(norm(x)).*cos(2*angle(x(:,2),x(:,1))/3+pi).*cos(angle(x(:,2),x(:,1))).*[norm(x)>.25 & norm(x)<.75];    
      fprintf(fid,'Singular function 4 ');
  end
  
                                 
  
  
  % Initialize mesh
  if M_flag == 1 
  Mesh.Coordinates = [0 0;1 0; 1 1; 0 1];
  Mesh.Elements = [1 2 3; 1 3 4];
  fprintf(fid,'in square domain\n');
  else
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  %Mesh = rotate(Mesh,pi/2);
  fprintf(fid,'in L-shaped domain\n');
  end
  
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);         
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = refine_REG(Mesh);
  
  fprintf(fid,'Mesh Width   L2 Norm of exact solution   L2 Difference        Order  \n');
  fprintf(fid,'==========   =========================   =================   =========\n');
  
  for i = 1:NREFS
  Mesh = refine_REG(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  
  % Assemble matrices and load vectors (discontinuous Lagrangian elements)
  
  QuadRule_1D = gauleg(0,1,2);
  QuadRule_2D = P3O3();
  SIGMA = @(P0,P1,varargin)alpha1/norm(P1-P0);  % Edge weight function
  [I1,J1,Avol_Curl] = assemMat_Vol_DG2(Mesh,@STIMA_Curl_LFE2);
  [I1,J1,Avol_Div] = assemMat_Vol_DG2(Mesh,@STIMA_Div_LFE2);
  
  [I2,J2,Jinn] = assemMat_Inn_DG2(Mesh,@STIMA_InnPen_DGLFE2,SIGMA);
  [I2,J2,Ainn_Curl] = assemMat_Inn_DG2(Mesh,@STIMA_Curl_Inn_DGLFE2,S);   
  [I2,J2,Ainn_Div] = assemMat_Inn_DG2(Mesh,@STIMA_Div_Inn_DGLFE2,S);
    
  [I3,J3,Jbnd] = assemMat_Bnd_DG2(Mesh,@STIMA_BndPen_Tangent_DGLFE2,SIGMA); 
  [I3,J3,Abnd_Curl] = assemMat_Bnd_DG2(Mesh,@STIMA_Curl_Bnd_DGLFE2,S);   
 
  Lvol = assemLoad_Vol_DG2(Mesh,@LOAD_Vol_DGLFE2,QuadRule_2D,G1,G2);
  Lbnd = assemLoad_Bnd_DG2(Mesh,@LOAD_Curl_Bnd_Tangent_DGLFE2,QuadRule_1D,S,SIGMA,UD1,UD2);
   
  % Create system matrix
       
  A = sparse([J1; J2; J3; J2; J3], ...
             [I1; I2; I3; I2; I3], ...
             [Avol_Curl+Avol_Div; Ainn_Curl+Ainn_Div; Abnd_Curl; Jinn; Jbnd]);
  L = Lvol+Lbnd;
  
  % Solve the linear system
  U1 = A\L;
  
  clear A L Avol_Curl Avol_Div Ainn_Curl Ainn_Div Abnd_Curl Jinn Jbnd I1 I2 I3 J1 J2 J3
  SIGMA = @(P0,P1,varargin)alpha2/norm(P1-P0);  % Edge weight function
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
  
  U2 = A\L;
  
  % Compute the errors
  L2Err(i) =  L2Err_DGLFE2(Mesh,U1-U2,P7O6(),O_Handle,O_Handle);
  Norm(i) = L2Err_DGLFE2(Mesh,zeros(size(U1,1),1),P7O6(),UD1,UD2);
  M_W(i) = get_MeshWidth(Mesh);
  if i==1
      fprintf(fid,'%8.4f   %20.4e       %17.4e          -\n',M_W(i),Norm(i),L2Err(i));
  else
      fprintf(fid,'%8.4f   %20.4e       %17.4e       %6.5f\n',M_W(i),Norm(i),L2Err(i),log(L2Err(i-1)/L2Err(i))/log(2));
  end

  clear A L U
  end
  fprintf(fid,'\n');
  end
  end
  fclose(fid);
  
  % Clear memory
    
  clear all
