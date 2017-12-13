% Run script for discontinuous Galerkin finite element solver
% Numerical results will be saved in the file WReg_mix.txt

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  fid = fopen('WReg_mix.txt','w');

  S = -1;        % Symmetric (+1) or antisymmetric (-1) discretization
  NREFS = 5;     % Number of red refinement steps
  U_Handle = @(x,varargin)ones(size(x,1),1);
  O_Handle = @(x,varargin)0*[ones(size(x,1),1) ones(size(x,1),1)];
  clear Mesh
  fprintf(fid,'Weak regularization version with edge element\n');  
  fprintf(fid,'\n');
  for M_flag = 1:2
  for F_flag = [1 3 4]
  clear Mesh 
  if F_flag == 1
      G1 = @(x,varargin)-2*ones(size(x,1),1);  
      G2 = @(x,varargin)2*ones(size(x,1),1);   
      UD1 = @(x,varargin)x(:,2).*(1+x(:,2));   
      UD2 = @(x,varargin)x(:,1).*(1-x(:,1));  
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
      fprintf(fid,'Singular function 1 ');
  elseif F_flag == 4
      G1 = @(x,varargin)0*ones(size(x,1),1);  
      G2 = @(x,varargin)0*ones(size(x,1),1);   
      UD1 = @(x,varargin)2/3*(x(:,1).^2+x(:,2).^2).^(-1/6).*sin(-angle(x(:,2),x(:,1))/3);    
      UD2 = @(x,varargin)2/3*(x(:,1).^2+x(:,2).^2).^(-1/6).*cos(-angle(x(:,2),x(:,1))/3);   
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
  F_Handle = @(x,varargin)[G1(x) G2(x)];
  GD_Handle = @(x,varargin)[UD1(x) UD2(x)];
                                 
  SIGMA = @(P0,P1,varargin)alpha/norm(P1-P0);  % Edge weight function
  
  % Initialize mesh
  if M_flag == 1 
  Mesh.Coordinates = [0 0;1 0; 1 1; 0 1];
  Mesh.Elements = [1 2 3; 1 3 4];
  fprintf(fid,'in square domain\n');
  else
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  fprintf(fid,'in L-shaped domain\n');
  end
  
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);         
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = refine_REG(Mesh);
  
  fprintf(fid,'Mesh Width   L2 Norm of exact solution   Relative L2 Error     Order  \n');
  fprintf(fid,'==========   =========================   =================   =========\n');
  
  for i = 1:NREFS
  Mesh = refine_REG(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  
  % Assemble matrices and load vectors (discontinuous Lagrangian elements)
  
  QuadRule_1D = gauleg(0,1,10);
  QuadRule_2D = P3O3();
  
        [IC,JC,C] = assemMat_W1F(Mesh,@STIMA_Curl_W1F,U_Handle,P7O6());
        B = assemMat_WRegW1F(Mesh,@STIMA_WReg_W1F);
        D = assemMat_LFE(Mesh,@MASS_Lump_LFE);
        
        nCoordinates = size(Mesh.Coordinates,1);
        Loc = get_BdEdges(Mesh);
        DEdges = Loc(Mesh.BdFlags(Loc) == -1);
        DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
        FreeDofs = setdiff(1:nCoordinates,DNodes);
        B = B(:,FreeDofs);
        D = D(FreeDofs,FreeDofs);
    
        T = B*inv(D)*transpose(B);
        [IT,JT,T] = find(T);
        A = sparse([IC;IT],[JC;JT],[C;T]);
        L = assemLoad_W1F(Mesh,P7O6(),F_Handle);

        % Incorporate Dirichlet boundary data

        [U,FreeDofs] = assemDir_W1F(Mesh,-1,GD_Handle,gauleg(0,1,1));
        L = L - A*U;

        % Solve the system

        U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        
  % Compute the errors
        L2Err(i) = L2Err_W1F(Mesh,U,P7O6(),GD_Handle);
        M_W(i) = get_MeshWidth(Mesh);
        Norm_Ex(i) = L2Err_W1F(Mesh,zeros(size(U,1),1),P7O6(),GD_Handle);
        Norm(i) = L2Err_W1F(Mesh,U,P7O6(),O_Handle);
       % L2Err(i) = L2Err(i)/Norm_Ex(i);
  if i==1
      fprintf(fid,'%8.4f & %20.4e   &   %17.4e    &     -   \\\\ \n',M_W(i),Norm(i),L2Err(i));
  else
      fprintf(fid,'%8.4f & %20.4e   &   %17.4e    &  %6.5f\\\\ \n',M_W(i),Norm(i),L2Err(i),log(L2Err(i-1)/L2Err(i))/log(2));
  end
  clear A L U
  end
  
    if M_flag==1 
      if F_flag==1
          save WReg11.mat M_W L2Err
      elseif F_flag==3
          save WReg13.mat M_W L2Err
      elseif F_flag==4
          save WReg14.mat M_W L2Err
      end
    else
      if F_flag==1
          save WReg21.mat M_W L2Err
      elseif F_flag==3
          save WReg23.mat M_W L2Err
      elseif F_flag==4
          save WReg24.mat M_W L2Err
      end
    end
  
  fprintf(fid,'\n');
  end
  end
  fclose(fid);
  % Clear memory
    
  clear all
