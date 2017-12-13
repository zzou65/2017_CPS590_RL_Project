% Run script for discontinuous Galerkin finite element solver
% Numerical results will be saved in the file SReg_LFE2.txt

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  fid = fopen('SReg_LFE2.txt','w');
  
  
  NREFS = 5;     % Number of red refinement steps
  U_Handle = @(x,varargin)ones(size(x,1),1);
  O_Handle = @(x,varargin)0*[ones(size(x,1),1) ones(size(x,1),1)];
  clear Mesh
  fprintf(fid,'Strong regularization version with nodal element\n');  
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
      UD1 = @(x,varargin)sing_fcn(x).*sin(-angle(x(:,2),x(:,1))/3);    
      UD2 = @(x,varargin)sing_fcn(x).*cos(-angle(x(:,2),x(:,1))/3);   
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
      UD1 = @(x,varargin)sing_fcn(x).*sin(angle(x(:,2),x(:,1))/3+pi).*[norm(x)<=.25]+sing_fcn(x).*sin(angle(x(:,2),x(:,1))/3+pi).*phi(norm(x)).*[norm(x)>.25 & norm(x)<.75]+norm(x).^(2/3).*phi1(norm(x)).*cos(2*angle(x(:,2),x(:,1))/3+pi).*sin(angle(x(:,2),x(:,1))).*[norm(x)>.25 & norm(x)<.75];    
      UD2 = @(x,varargin)-sing_fcn(x).*cos(angle(x(:,2),x(:,1))/3+pi).*[norm(x)<=.25]-sing_fcn(x).*cos(angle(x(:,2),x(:,1))/3+pi).*phi(norm(x)).*[norm(x)>.25 & norm(x)<.75]-norm(x).^(2/3).*phi1(norm(x)).*cos(2*angle(x(:,2),x(:,1))/3+pi).*cos(angle(x(:,2),x(:,1))).*[norm(x)>.25 & norm(x)<.75];    
      fprintf(fid,'Singular function 4 ');
  end
  F_Handle = @(x,varargin)[G1(x) G2(x)];
  GD_Handle = @(x,varargin)[UD1(x) UD2(x)];
                                 
  SIGMA = @(P0,P1,varargin)alpha/norm(P1-P0);  % Edge weight function
  
  % Initialize mesh
  if M_flag == 1 
  Mesh.Coordinates = [0 0;1 0; 1 1; 0 1];
  Mesh.Elements = [1 2 3; 1 3 4];
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);         
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  

  fprintf(fid,'in square domain\n');
  else
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  Mesh = add_Edges(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Loc = get_BdEdges(Mesh);
  BdFlags = [-1 -2 -3 -4 -5 -6];
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc)= -1;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  fprintf(fid,'in L-shaped domain\n');
  end
  
  Mesh = refine_REG(Mesh);
  
  fprintf(fid,'Mesh Width   L2 Norm of exact solution   Relative L2 Error     Order  \n');
  fprintf(fid,'==========   =========================   =================   =========\n');
  
  for i = 1:NREFS
  Mesh = refine_REG(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  
  % Assemble matrices and load vectors (discontinuous Lagrangian elements)
  
  QuadRule_1D = gauleg(0,1,2);
  QuadRule_2D = P3O3();
  
        nCoordinates = size(Mesh.Coordinates,1);
        [IC,JC,C] = assemMat_LFE2(Mesh,@STIMA_Curl_LFE2,U_Handle,P7O6());
        if M_flag == 1 | M_flag == 2
            [ID,JD,D] = assemMat_LFE2(Mesh,@STIMA_Div_LFE2,U_Handle,P7O6());
        else 
            [ID,JD,D] = assemMat_LFE2(Mesh,@STIMA_Div_Weighted2_LFE2,U_Handle,P7O6());
            D = 2*D;
        end
        [U,g,FreeDofs,IB,JB,B] = assemDir_StrRegLFE2(Mesh,-1,GD_Handle);
        A = sparse([IC;ID;IB+2*nCoordinates;JB],[JC;JD;JB;IB+2*nCoordinates],[C;D;B;B]);
        l = assemLoad_LFE2(Mesh,P7O6(),F_Handle);
        L = [l;g];
        L = L - A*U;

        % Solve the system

        U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        U = U(1:2*nCoordinates);

      %  if i ==NREFS
      %        figure
      %  plot_LFE2(U(1:end/2),Mesh);
      %  figure
      %  plot_LFE2(U(end/2+1:end),Mesh);
      %  end
        
        % Compute the errors
        L2Err(i) = L2Err_LFE2(Mesh,U,P7O6(),GD_Handle);
        M_W(i) = get_MeshWidth(Mesh);
        Norm_Ex(i) = L2Err_LFE2(Mesh,zeros(size(U,1),1),P7O6(),GD_Handle);
        Norm(i) = L2Err_LFE2(Mesh,U,P7O6(),O_Handle);
        %L2Err(i) = L2Err(i)/Norm_Ex(i);
  if i==1
      fprintf(fid,'%8.4f & %20.4e   &   %17.4e    &     -   \\\\ \n',M_W(i),Norm(i),L2Err(i));
  else
      fprintf(fid,'%8.4f & %20.4e   &   %17.4e    &  %6.5f\\\\ \n',M_W(i),Norm(i),L2Err(i),log(L2Err(i-1)/L2Err(i))/log(2));
  end
save SReg_Sol.mat Mesh U;
  clear A L U
  end
  
    if M_flag==1 
      if F_flag==1
          save SReg11.mat M_W L2Err
      elseif F_flag==3
          save SReg13.mat M_W L2Err
      elseif F_flag==4
          save SReg14.mat M_W L2Err
      end
    else
      if F_flag==1
          save SReg21.mat M_W L2Err
      elseif F_flag==3
          save SReg23.mat M_W L2Err
      elseif F_flag==4
          save SReg24.mat M_W L2Err
      end
    end
    
  fprintf(fid,'\n');
  end
  end
  fclose(fid);

  % Clear memory
    
  clear all
