function [] = main6()
% convergence of local multigrid on variable refined meshes

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize mesh

  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  
  % apply local multigrid
  
  maxit = 50;
  tol = 0.015;
  m = [1 1];
  [U,Mesh,flag,iter,ndof,err,timer,dofslvl,A,b,FDofs] = ...
    locmg_solve([],Mesh,@f_LShap,@g_D_LShap,'h',0.5,1,1,m,@gs_smooth,tol,maxit);

  % calculate exact solution and check convergence
  
  U_ex = U;
  U_ex(FDofs) = 0;
  U_ex(FDofs) = A(FDofs,FDofs)\(b(FDofs)-A(FDofs,:)*U_ex);
  err = sqrt((U-U_ex)'*A*(U-U_ex));
  disc_err = H1SErr_LFE(Mesh,U_ex,P7O6(),@grad_uex_LShap);
  fprintf('disc. error : %g    trunc. error : %g\n',disc_err,err);
  
  % calculate number of 1-dim ssc at each iteration
  
  numssc = zeros(1,iter);
  for i=1:iter
    numssc(i:iter) = numssc(i:iter) + (m(1)+m(end))*sum(dofslvl(:,i));
  end
  
  % crop data
  
  iter1 = ceil(0.3*iter);
  numssc = numssc(iter1:iter);
  ndof = ndof(iter1:iter);
  
  % plot number of subspace corrections vs number of dofs
  
  figure;
  plot(ndof,numssc,'o');
  set(gca,'XScale','log','YScale','log');
  grid('on');
  hold on;
  
  title('{\bf Complexity of Local Mutligrid}');
  xlabel('{\bf degrees of freedom}');
  ylabel('{\bf one-dimensonal subspace corrections}')
  
  p = polyfit(log(ndof),log(numssc),1);
  x = get(gca,'XLim');
  y = get(gca,'YLim');
  plot(x,exp(polyval(p,log(x))));
  set(gca,'YLim',y);
  add_Slope(gca,'NorthWest',p(1));
  
return