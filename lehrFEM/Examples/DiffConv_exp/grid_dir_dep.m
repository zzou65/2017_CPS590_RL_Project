% Initialize constants
  
 clear Mesh;
 
  JIG=1;
  d1=2;                    %impact of supg modification
  d2=2;
  a=10^-10                %amount of diffusivity
  NREFS =3;
  
  % solution u
  d.U_EX_Handle=@(x,d,a)(x(:,1)<=(ones(size(x,1),1)*0.5));
  % solution grad u
  d.GRAD_U_EX_Handle=...
            @(x,varargin)0;                    
  % -a*div grad u+b*v grad u
  d.SOL_Handle=@(x,elemflag,a,b)zeros(size(x,1),1);
  d.Coordinates=[0 0; 1 0; 1 1; 0 1];
  d.Elements=[1 2 4;2 3 4];
  d.boundtype=[-1 -1 -1 -1]; 
  
  nTheta=5;
  Thetas=linspace(pi/16,15*pi/16,nTheta);
%    nTheta=1;
%    Thetas=pi/2;
  QuadRule = P7O6(); 
  
  Mesh.Coordinates =d.Coordinates; 
  Mesh.Elements = d.Elements;
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = d.boundtype;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  err=zeros(NREFS,5);
  err1=zeros(NREFS,5);
  h=zeros(NREFS,1);
  Dofs=zeros(NREFS,1);
  plot_Mesh(Mesh,'a');
  print -depsc ~/griddep_mesh.eps
  for i = 1:NREFS
   
   %refine Mesh   
   Mesh = refine_REG(Mesh);
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
  end       
   
  %Laplace
  A = assemMat_QFE(New_Mesh,@STIMA_LAPL_QFE_quad,P3O3());
  L = assemLoad_QFE(New_Mesh,P3O3(),d.SOL_Handle,a);
      
  % theta dependend terms
  
  %xi=(0:0.01*get_MeshWidth(New_Mesh):sqrt(2))-sqrt(2)/2;
  xi=(0:0.001:sqrt(2))-sqrt(2)/2;
  np=size(xi,2);
  U_1=zeros(np,nTheta); U_2=zeros(np,nTheta);
  U_3=zeros(np,nTheta); U_4=zeros(np,nTheta);
  U_5=zeros(np,nTheta); U_6=zeros(np,nTheta);
  
%   Thetas=pi/2;
  for i=1:nTheta 
      
   theta=Thetas(i);
   % velocity v pi/4<=theta <=3/4pi 
   d.V_Handle=...
            @(x,varargin)[cos(theta)*ones(size(x,1),1) sin(theta)*ones(size(x,1),1)];
   L_supg=assemLoad_QFE_SUPG(New_Mesh,P3O3(),d.V_Handle,d.SOL_Handle,a,d1,d2);
   
   weights1=[0 0 0 1/3 1/3 1/3 0];
   weights2=[1/12 1/12 1/12 0 0 0 3/4];
   weights3=[0.05 0.05 0.05 4/30 4/30 4/30, 0.45];
   weights4=[1/3 1/3 1/3 0 0 0 0];
   [M1,UP1,C1]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights1);
   [M2,UP2,C2]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights2);
   [M3,UP3,C3]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights3);
   [M4,UP4,C4]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights4);
   B_supg=assemMat_QFE(New_Mesh,@STIMA_SUPG_QFE,P3O3(), d.V_Handle,a,d1,d2);
   B=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,P3O3());
   A_u1=a*A+(M1*UP1);
   A_u2=a*A+(M2*UP2+C2);
   A_u3=a*A+(M3*UP3+C3);
   A_u4=a*A+(M4*UP4+C4);
   A_supg=a*A+B+B_supg;
   A_s=a*A+B;
   
   %Incorporate Dirichlet boundary data
   
   [U_s,FreeDofs] = assemDir_QFE(New_Mesh,[-1 -2],d.U_EX_Handle,a);
   U_u1=U_s;
   U_u2=U_s;
   U_u3=U_s;
   U_u4=U_s;
   U_supg=U_s;
    
   L_u1 = L - A_u1*U_u1;
   L_u2 = L - A_u2*U_u2;
   L_u3 = L - A_u3*U_u3;
   L_u4 = L - A_u4*U_u4;
   L_s = L - A_s*U_s;
   L_supg = L+L_supg - A_supg*U_supg;
   
   % Solve the linear system
 
   U_u1(FreeDofs) = A_u1(FreeDofs,FreeDofs)\L_u1(FreeDofs);
   U_u2(FreeDofs) = A_u2(FreeDofs,FreeDofs)\L_u2(FreeDofs);
   U_u3(FreeDofs) = A_u3(FreeDofs,FreeDofs)\L_u3(FreeDofs);
   U_u4(FreeDofs) = A_u4(FreeDofs,FreeDofs)\L_u4(FreeDofs);
   U_s(FreeDofs) = A_s(FreeDofs,FreeDofs)\L_s(FreeDofs);
   U_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);
   %
   xstart=[1/2 1/2]-1/(2*sin(theta))*[sin(theta), -cos(theta)];
   xend=[1/2 1/2]+1/(2*sin(theta))*[sin(theta), -cos(theta)];

   [x,y]=plotLine_QFE([U_u4 U_u1 U_u2 U_u3 U_s U_supg],New_Mesh,xstart,xend); 
   yi=interp1(x-1/(2*sin(theta)),y,xi);
   U_1(1:np,i)=yi(1:np,1);
   U_2(1:np,i)=yi(1:np,2);
   U_3(1:np,i)=yi(1:np,3);
   U_4(1:np,i)=yi(1:np,4);
   U_5(1:np,i)=yi(1:np,5);
   U_6(1:np,i)=yi(1:np,6);
  % plot_QFE(U_supg,New_Mesh);
  
  I_3=inv(full(A_u3(FreeDofs,FreeDofs)));
  norm(I_3,inf)
  end
  get_MeshWidth(New_Mesh)
  [X,Y]=meshgrid(xi,Thetas);
  XLim=[min(xi)-0.05 max(xi)+0.05];
  Ylim=[min(Thetas)-0.05 max(Thetas)+0.05];
  figure; plot3(X',Y',U_1,'LineWidth',3); grid on; 
  set(gca,'XLim',XLim,'YLim',YLim,'ZLim',[min(min(U_1))-0.05 max(max(U_1))+0.05],'FontSize',13);
  xlabel('{\bf profile line}'); ylabel('{\bf \Theta}'); zlabel('{\bf profile}');
  print -depsc ~/griddep_1.eps
  figure; plot3(X',Y',U_2,'LineWidth',3); grid on;
  set(gca,'XLim',XLim,'YLim',YLim,'ZLim',[min(min(U_2))-0.05 max(max(U_2))]+0.05,'FontSize',13);
  xlabel('{\bf profile line}'); ylabel('{\bf \Theta}'); zlabel('{\bf profile}');
  print -depsc ~/griddep_2.eps
  figure; plot3(X',Y',U_3,'LineWidth',3); grid on;
  set(gca,'XLim',XLim,'YLim',YLim,'ZLim',[min(min(U_3))-0.05 max(max(U_3))]+0.05,'FontSize',13);
  xlabel('{\bf profile line}'); ylabel('{\bf \Theta}'); zlabel('{\bf profile}');
  print -depsc ~/griddep_3.eps
  figure; plot3(X',Y',U_4,'LineWidth',3); grid on;
  set(gca,'XLim',XLim,'YLim',YLim,'ZLim',[min(min(U_4))-0.05 max(max(U_4))]+0.05,'FontSize',13);
  xlabel('{\bf profile line}'); ylabel('{\bf \Theta}'); zlabel('{\bf profile}');
  print -depsc ~/griddep_4.eps
  figure; plot3(X',Y',U_5,'LineWidth',3); grid on;
  set(gca,'XLim',XLim,'YLim',YLim,'ZLim',[min(min(U_5))-0.05 max(max(U_5))]+0.05,'FontSize',13);  
  xlabel('{\bf profile line}'); ylabel('{\bf \Theta}'); zlabel('{\bf profile}');
  print -depsc ~/griddep_s.eps
  figure; plot3(X',Y',U_6,'LineWidth',3); grid on;
  set(gca,'XLim',XLim,'YLim',YLim,'ZLim',[min(min(U_6))-0.05 max(max(U_6))]+0.05,'FontSize',13);
  xlabel('{\bf profile line}'); ylabel('{\bf \Theta}'); zlabel('{\bf profile}');
  print -depsc ~/griddep_supg.eps
  figure; 
  subplot(2,3,1); plot(xi,U_1);
  set(gca,'XLim',XLim,'YLim',[min(min(U_1))-0.05 max(max(U_1))+0.05]);
  xlabel('{\bf profile line}'); ylabel('{\bf profile}');
  subplot(2,3,2); plot(xi,U_2);
  set(gca,'XLim',XLim,'YLim',[min(min(U_2))-0.05 max(max(U_2))+0.05]);
  xlabel('{\bf profile line}'); ylabel('{\bf profile}');
  subplot(2,3,3); plot(xi,U_3);
  set(gca,'XLim',XLim,'YLim',[min(min(U_3))-0.05 max(max(U_3))+0.05]);
  xlabel('{\bf profile line}'); ylabel('{\bf profile}');
  subplot(2,3,4); plot(xi,U_4);
  set(gca,'XLim',XLim,'YLim',[min(min(U_4))-0.05 max(max(U_4))]+0.05);
  xlabel('{\bf profile line}'); ylabel('{\bf profile}');
  subplot(2,3,5); plot(xi,U_5);
  set(gca,'XLim',XLim,'YLim',[min(min(U_5))-0.05 max(max(U_5))]+0.05);
  xlabel('{\bf profile line}'); ylabel('{\bf profile}');
  subplot(2,3,6); plot(xi,U_6);
  set(gca,'XLim',XLim,'YLim',[min(min(U_6))-0.05 max(max(U_6))]+0.05);
  xlabel('{\bf profile line}'); ylabel('{\bf profile}');
 
clear all;