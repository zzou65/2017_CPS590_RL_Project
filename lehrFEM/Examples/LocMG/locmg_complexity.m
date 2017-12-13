function [] = locmg_complexity()
% complexity of local (adaptive) multigrid solver
%
%   This code plots the time required for various parts of the locmg_solve
%   code (local multigrid) against the number of degrees of freedom.

%   Copyright 2006-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % Initialize constants

  F_HANDLE = @f_LShap;
  GD_HANDLE = @g_D_LShap;
  CMesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  U0 = zeros(size(CMesh.Coordinates,1),1);
  M = [2,1];
  CYC = 1;
  SMOOTHER = @gs_smooth;
  TOL = (1e-2)*[8,6,4,3,2,1.6];
  MAXIT = 40;
  THETA = 0.5;
  ERREST = 'c';
  ITLVL = 1;

  % Run multigrid solver

  ndofs = zeros(size(TOL));
  time_total = zeros(size(TOL)); 
  time_mg = zeros(size(TOL)); 
  time_errest = zeros(size(TOL)); 
  time_mesh = zeros(size(TOL)); 
  time_stima = zeros(size(TOL));
  time_ref = zeros(size(TOL));

  for k=1:length(TOL)

    % calculate number of degrees of freedom

    [U,Mesh,flag,iter,ndof,err,timer] = locmg_solve(U0,CMesh,F_HANDLE,GD_HANDLE,ERREST,THETA,ITLVL,CYC,M,SMOOTHER,TOL(k),MAXIT);
    time = timer.total;
    ndofs(k) = ndof(iter);

    % time method accordingly

    num = max(floor(10/time),1);
    temp = repmat(timer,1,num);
    t = cputime;
    for i=1:num
      [U,Mesh,flag,iter,ndof,err,temp(i)] = locmg_solve(U0,CMesh,F_HANDLE,GD_HANDLE,ERREST,THETA,ITLVL,CYC,M,SMOOTHER,TOL(k),MAXIT);
    end
    time_total(k) = (cputime-t)/num;
    tempSum = temp(1);
    for i=2:num
      tempSum.total = tempSum.total + temp(i).total;
      tempSum.mg = tempSum.mg + temp(i).mg;
      tempSum.errest = tempSum.errest + temp(i).errest;
      tempSum.mesh = tempSum.mesh + temp(i).mesh;
      tempSum.stima = tempSum.stima + temp(i).stima;
      tempSum.ref = tempSum.ref + temp(i).ref;
    end
    time_mg(k) = tempSum.mg/num;
    time_errest(k) = tempSum.errest/num;
    time_mesh(k) = tempSum.mesh/num;
    time_stima(k) = tempSum.stima/num;
    time_ref(k) = tempSum.ref/num;
  end

  % Calculate complexities

  p_total = polyfit(log(ndofs),log(time_total),1);
  p_mg = polyfit(log(ndofs),log(time_mg),1);
  p_errest = polyfit(log(ndofs),log(time_errest),1);
  p_mesh = polyfit(log(ndofs),log(time_mesh),1);
  p_stima = polyfit(log(ndofs),log(time_stima),1);
  p_ref = polyfit(log(ndofs),log(time_ref),1);

  % Plot time against number of degrees of freedom

  figure;
  plot(ndofs,time_total,'bo');
  hold on;
  plot(ndofs,time_mg,'go');
  plot(ndofs,time_errest,'co');
  plot(ndofs,time_mesh,'mo');
  plot(ndofs,time_stima,'ro');
  plot(ndofs,time_ref,'ko');
  set(gca,'XScale','log','YScale','log');
  grid('on');

  title('{\bf complexity of local multigrid implementation}');
  xlabel('{\bf degrees of freedom}');
  ylabel('{\bf time}')

  legend(sprintf('total (%3.2f)',p_total(1)),...
         sprintf('MG cycles (%3.2f)',p_mg(1)),...
         sprintf('error est. (%3.2f)',p_errest(1)),...
         sprintf('mesh (%3.2f)',p_mesh(1)),...
         sprintf('assembly (%3.2f)',p_stima(1)),...
         sprintf('refinement (%3.2f)',p_ref(1)),...
         'Location','NorthWest');

  x = get(gca,'XLim');
  y = get(gca,'YLim');

  plot(x,exp(polyval(p_total,log(x))),'b');
  plot(x,exp(polyval(p_mg,log(x))),'g');
  plot(x,exp(polyval(p_errest,log(x))),'c');
  plot(x,exp(polyval(p_mesh,log(x))),'m');
  plot(x,exp(polyval(p_stima,log(x))),'r');
  plot(x,exp(polyval(p_ref,log(x))),'k');

  set(gca,'YLim',y);

return