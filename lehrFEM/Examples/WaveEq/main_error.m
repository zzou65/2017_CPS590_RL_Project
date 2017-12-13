% 1D Wave equation with newmark anf theta scheme 
%
% Copyright 2008-2008 Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  A=0;                                                     % interval beginning 
  L = 1;                                                    % Length of the interval
  T = 1;                                                    % Final time
  U1 = @(x,varargin)x.*(1-x);                     % Initial data
  U2 = @(x,varargin)sin(pi*x);                     % Initial data
  N=8;
  NPTS = 2.^((1:N)+1);                                                % Number of points
  NSTEPS = 2.^((1:N)+1);                                             % Number of time steps 
 
  % Schemes are identical if gamma=1/2 beta=1/4 and theta=1/2;
  gamma=1/2; 
  beta=1;%1/4;
  theta=1/2;
  
  err=zeros(N,2);
  tic
  % error calculation
  for i=1:N
    err(i,1)=theta_error(A,L,T,U1,NPTS(i),NSTEPS(i),theta);
    err(i,2)=theta_error(A,L,T,U2,NPTS(i),NSTEPS(i),theta);  
    err(i,3)=newmark_error(A,L,T,U1,NPTS(i),NSTEPS(i),beta,gamma);
    err(i,4)=newmark_error(A,L,T,U2,NPTS(i),NSTEPS(i),beta,gamma);
  end
  toc
 
 xstep=1./(NPTS-1);
 errL2=err;
 figure;
 plot(xstep,errL2(:,1),'b',xstep,errL2(:,2),'r',xstep,errL2(:,3),'g',xstep,errL2(:,4),'c');
 set(gca,'XScale','log','YScale','log'); 
 xlabel('\bf Meshwidth');
 ylabel('\bf L^2-error');
 grid('on');
 legend('x-x^2, \theta-scheme','sin(pi*x), \theta-scheme','x-x^2, Newark scheme','sin(pi*x), Newark scheme');
 p = polyfit(log(xstep(end-3:end)),log(errL2(end-3:end,1)'),1);
 add_Slope(gca,'North',p(1),'b');
 p = polyfit(log(xstep(end-3:end)),log(errL2(end-3:end,2)'),1);
 add_Slope(gca,'NorthWest',p(1),'r');
 p = polyfit(log(xstep(end-3:end)),log(errL2(end-3:end,3)'),1);
 add_Slope(gca,'South',p(1),'g');
 p = polyfit(log(xstep(end-3:end)),log(errL2(end-3:end,4)'),1);
 add_Slope(gca,'SouthWest',p(1),'c');
 
 print -depsc plot.eps
  
