 % Initialize constants
  clear Mesh
 
  d1=0.1;                    %impact of supg modification
  d2=0.1;
  a=1                %amount of diffusivity
  NREFS =1;
  V_Handle=@(x,varargin)ones(size(x,1),1)*[cos(pi/6) sin(pi/6)];
  QuadRule = P7O6(); 
  
  Mesh.Coordinates =[-1 -1; 1 -1; 1 1; -1 1];
  Mesh.Elements = [1 2 4; 2 3 4];
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  
  %refine Mesh   
  Mesh = refine_REG(Mesh);
  Mesh=add_Edge2Elem(Mesh);
        
  New_Mesh = Mesh;      
        
   %Laplace
   A = assemMat_QFE(New_Mesh,@STIMA_Lapl_QFE);
   weights1=[0 0 0 1/3 1/3 1/3 0];
   weights2=[1/12 1/12 1/12 0 0 0 3/4];
   weights3=[0.05 0.05 0.05 4/30 4/30 4/30, 0.45];
   weights4=[1/3 1/3 1/3 0 0 0 0];
   [M1,UP1,C1]=assemMat_UpQFE(New_Mesh,V_Handle,weights1);
   [M2,UP2,C2]=assemMat_UpQFE(New_Mesh,V_Handle,weights2);
   [M3,UP3,C3]=assemMat_UpQFE(New_Mesh,V_Handle,weights3);
   [M4,UP4,C4]=assemMat_UpQFE(New_Mesh,V_Handle,weights4);
   B_supg=assemMat_QFE(New_Mesh,@STIMA_SUPG_QFE,P7O6(), V_Handle,a,d1,d2);
   B=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,V_Handle,P7O6());
   A_u1=M1*UP1;
   A_u2=M2*UP2+C2;
   A_u3=M3*UP3+C3;
   A_u4=M4*UP4+C4;
   A_supg=B+B_supg;
   nodemap=[4 19 9 15 3 17 18 24 25 16 6 22 7 23 8 11 20 21 12 14 1 10 5 13 2];
   
   % Laplace stencils
   Lapl_1=reshape(full(A(7,nodemap)),5,5)';
   Lapl_2=reshape(full(A(24,nodemap)),5,5)';
   Lapl_3=reshape(full(A(23,nodemap)),5,5)';
   Lapl_4=reshape(full(A(25,nodemap)),5,5)';
   
   %Convection stencils
   Conv_1=reshape(full(B(7,nodemap)),5,5)';
   Conv_2=reshape(full(B(24,nodemap)),5,5)';
   Conv_3=reshape(full(B(23,nodemap)),5,5)';
   Conv_4=reshape(full(B(25,nodemap)),5,5)';
   
   %upwind stencils 
   %1) using vertices
   Up1_1=reshape(full(A_u4(7,nodemap)),5,5)';
   Up1_2=reshape(full(A_u4(24,nodemap)),5,5)';
   Up1_3=reshape(full(A_u4(23,nodemap)),5,5)';
   Up1_4=reshape(full(A_u4(25,nodemap)),5,5)';
   
   %2) using midpoints
   Up2_1=reshape(full(A_u1(7,nodemap)),5,5)';
   Up2_2=reshape(full(A_u1(24,nodemap)),5,5)';
   Up2_3=reshape(full(A_u1(23,nodemap)),5,5)';
   Up2_4=reshape(full(A_u1(25,nodemap)),5,5)';
   
   %3) vertices and barycenter
   Up3_1=reshape(full(A_u2(7,nodemap)),5,5)';
   Up3_2=reshape(full(A_u2(24,nodemap)),5,5)';
   Up3_3=reshape(full(A_u2(23,nodemap)),5,5)';
   Up3_4=reshape(full(A_u2(25,nodemap)),5,5)';
   
   %4) vertices, midpoints and barycenter
   Up4_1=reshape(full(A_u3(7,nodemap)),5,5)';
   Up4_2=reshape(full(A_u3(24,nodemap)),5,5)';
   Up4_3=reshape(full(A_u3(23,nodemap)),5,5)';
   Up4_4=reshape(full(A_u3(25,nodemap)),5,5)';
   
   %supg-stabilizaton
   SUPG_1=reshape(full(B_supg(7,nodemap)),5,5)';
   SUPG_2=reshape(full(B_supg(24,nodemap)),5,5)';
   SUPG_3=reshape(full(B_supg(23,nodemap)),5,5)';
   SUPG_4=reshape(full(B_supg(25,nodemap)),5,5)';
   
   k1min=-2; k1step=0.01; k1max=2;
   k2min=-2; k2step=0.01; k2max=2;
   p=3    % peclet  number  
   [k1, k2]=meshgrid(k1min:k1step:k1max , k2min:k2step:k2max);
   
   % standard
   stencil1=Lapl_1+p*Conv_1;
   stencil2=Lapl_2+p*Conv_2;
   stencil3=Lapl_3+p*Conv_3;
   stencil4=Lapl_4+p*Conv_4;
   e1=zeros(size(k1)); e2=zeros(size(k1));
   e2=zeros(size(k1)); e4=zeros(size(k1)); 
   n=size(k1,1);
   for i=1:1:n
       for j=1:1:n
            S=symbolQFE(stencil1, stencil2, stencil3, stencil4, p, k1(i,j), k2(i,j));
            e=eig(S);
            e1(i,j)=e(1); e2(i,j)=e(2); e3(i,j)=e(3); e4(i,j)=e(4);
       end
   end
   figure('Name','standard scheme','NumberTitle','off');
   subplot(2,2,1);  [C,h]=contour(k1,k2,real(e1)); set(h,'ShowText','on');
   subplot(2,2,2);  [C,h]=contour(k1,k2,real(e2)); set(h,'ShowText','on');
   subplot(2,2,3);  [C,h]=contour(k1,k2,real(e3)); set(h,'ShowText','on');
   subplot(2,2,4);  [C,h]=contour(k1,k2,real(e4)); set(h,'ShowText','on');
   
   % upwind using vertices  
   stencil1=Lapl_1+p*Up1_1;
   stencil2=Lapl_2+p*Up1_2;
   stencil3=Lapl_3+p*Up1_3;
   stencil4=Lapl_4+p*Up1_4;
   e1=zeros(size(k1)); e2=zeros(size(k1));
   e2=zeros(size(k1)); e4=zeros(size(k1)); 
   n=size(k1,1);
   for i=1:1:n
       for j=1:1:n
            S=symbolQFE(stencil1, stencil2, stencil3, stencil4, p, k1(i,j), k2(i,j));
            e=eig(S);
            e1(i,j)=e(1); e2(i,j)=e(2); e3(i,j)=e(3); e4(i,j)=e(4);
       end
   end
   figure('Name','upwind vertices','NumberTitle','off');  
   subplot(2,2,1);  [C,h]=contour(k1,k2,real(e1)); set(h,'ShowText','on');
   subplot(2,2,2);  [C,h]=contour(k1,k2,real(e2)); set(h,'ShowText','on');
   subplot(2,2,3);  [C,h]=contour(k1,k2,real(e3)); set(h,'ShowText','on');
   subplot(2,2,4);  [C,h]=contour(k1,k2,real(e4)); set(h,'ShowText','on');
   
   % upwind using midpoints  
   stencil1=Lapl_1+p*Up2_1;
   stencil2=Lapl_2+p*Up2_2;
   stencil3=Lapl_3+p*Up2_3;
   stencil4=Lapl_4+p*Up2_4;
   e1=zeros(size(k1)); e2=zeros(size(k1));
   e2=zeros(size(k1)); e4=zeros(size(k1)); 
   n=size(k1,1);
   for i=1:1:n
       for j=1:1:n
            S=symbolQFE(stencil1, stencil2, stencil3, stencil4, p, k1(i,j), k2(i,j));
            e=eig(S);
            e1(i,j)=e(1); e2(i,j)=e(2); e3(i,j)=e(3); e4(i,j)=e(4);
       end
   end
   figure('Name','upwind using  midpoints ','NumberTitle','off');  
   subplot(2,2,1);  [C,h]=contour(k1,k2,real(e1)); set(h,'ShowText','on');
   subplot(2,2,2);  [C,h]=contour(k1,k2,real(e2)); set(h,'ShowText','on');
   subplot(2,2,3);  [C,h]=contour(k1,k2,real(e3)); set(h,'ShowText','on');
   subplot(2,2,4);  [C,h]=contour(k1,k2,real(e4)); set(h,'ShowText','on');
   
   % upwind using vertices, barycenter  
   stencil1=Lapl_1+p*Up3_1;
   stencil2=Lapl_2+p*Up3_2;
   stencil3=Lapl_3+p*Up3_3;
   stencil4=Lapl_4+p*Up3_4;
   e1=zeros(size(k1)); e2=zeros(size(k1));
   e2=zeros(size(k1)); e4=zeros(size(k1)); 
   n=size(k1,1);
   for i=1:1:n
       for j=1:1:n
            S=symbolQFE(stencil1, stencil2, stencil3, stencil4, p, k1(i,j), k2(i,j));
            e=eig(S);
            e1(i,j)=e(1); e2(i,j)=e(2); e3(i,j)=e(3); e4(i,j)=e(4);
       end
   end
   figure('Name','upwind using vertices barycenter','NumberTitle','off');  
   subplot(2,2,1);  [C,h]=contour(k1,k2,real(e1)); set(h,'ShowText','on');
   subplot(2,2,2);  [C,h]=contour(k1,k2,real(e2)); set(h,'ShowText','on');
   subplot(2,2,3);  [C,h]=contour(k1,k2,real(e3)); set(h,'ShowText','on');
   subplot(2,2,4);  [C,h]=contour(k1,k2,real(e4)); set(h,'ShowText','on');
   
   % upwind using vertices, midpoints, barycenter  
   stencil1=Lapl_1+p*Up4_1;
   stencil2=Lapl_2+p*Up4_2;
   stencil3=Lapl_3+p*Up4_3;
   stencil4=Lapl_4+p*Up4_4;
   e1=zeros(size(k1)); e2=zeros(size(k1));
   e2=zeros(size(k1)); e4=zeros(size(k1)); 
   n=size(k1,1);
   for i=1:1:n
       for j=1:1:n
            S=symbolQFE(stencil1, stencil2, stencil3, stencil4, p, k1(i,j), k2(i,j));
            e=eig(S);
            e1(i,j)=e(1); e2(i,j)=e(2); e3(i,j)=e(3); e4(i,j)=e(4);
       end
   end
   figure('Name','upwind using vertices, midpoints, barycenter','NumberTitle','off');  
   subplot(2,2,1);  [C,h]=contour(k1,k2,real(e1)); set(h,'ShowText','on');
   subplot(2,2,2);  [C,h]=contour(k1,k2,real(e2)); set(h,'ShowText','on');
   subplot(2,2,3);  [C,h]=contour(k1,k2,real(e3)); set(h,'ShowText','on');
   subplot(2,2,4);  [C,h]=contour(k1,k2,real(e4)); set(h,'ShowText','on');
   
   % SUPG
   stencil1=Lapl_1+p*Conv_1+p*SUPG_1;
   stencil2=Lapl_2+p*Conv_2+p*SUPG_2;
   stencil3=Lapl_3+p*Conv_3+p*SUPG_3;
   stencil4=Lapl_4+p*Conv_4+p*SUPG_4;
   e1=zeros(size(k1)); e2=zeros(size(k1));
   e2=zeros(size(k1)); e4=zeros(size(k1)); 
   n=size(k1,1);
   for i=1:1:n
       for j=1:1:n
            S=symbolQFE(stencil1, stencil2, stencil3, stencil4, p, k1(i,j), k2(i,j));
            e=eig(S);
            e1(i,j)=e(1); e2(i,j)=e(2); e3(i,j)=e(3); e4(i,j)=e(4);
       end
   end
   figure('Name','SUPG scheme','NumberTitle','off');
   subplot(2,2,1);  [C,h]=contour(k1,k2,real(e1)); set(h,'ShowText','on');
   subplot(2,2,2);  [C,h]=contour(k1,k2,real(e2)); set(h,'ShowText','on');
   subplot(2,2,3);  [C,h]=contour(k1,k2,real(e3)); set(h,'ShowText','on');
   subplot(2,2,4);  [C,h]=contour(k1,k2,real(e4)); set(h,'ShowText','on');
   clear all;
