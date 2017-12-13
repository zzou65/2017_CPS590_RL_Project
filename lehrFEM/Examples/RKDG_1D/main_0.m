% Test routine for Legendre polynomials.
  
% Copyright 2007-2007 Patrick Meury & Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize grid points
  
  x = transpose(-1:0.01:1);
  
  % Evaluate Legendre polynomials (analytic values)
  
  P_0 = @(t)ones(size(t));
  P_1 = @(t)t;
  P_2 = @(t)(3*t.^2 - 1)/2;
  P_3 = @(t)(5*t.^3 - 3*t)/2;
  P_4 = @(t)(35*t.^4 - 30*t.^2 + 3)/8;
  
  DP_0 = @(t)zeros(size(t));
  DP_1 = @(t)ones(size(t));
  DP_2 = @(t)3*t;
  DP_3 = @(t)(15*t.^2 - 3)/2;
  DP_4 = @(t)(140*t.^3 - 60*t)/8;
  
  % Evaluate Legendre polynomials (recurrence relation)
  
  L = shap_Leg_1D(x,4)
  DL = grad_shap_Leg_1D(x,4)
  
  % Generate plots
  
  fig = figure('Name','Legendre polynomials (analytic values)');
  plot(x,P_0(x),'r-', ...
       x,P_1(x),'b-', ...
       x,P_2(x),'g-', ...
       x,P_3(x),'c-', ...
       x,P_4(x),'m-');
  title('{\bf Legendre polynomials (analytic values)}');
  
  fig = figure('Name','Legendre polynomials (recurrence relation)');
  plot(x,L(:,1),'r-', ...
       x,L(:,2),'b-', ...
       x,L(:,3),'g-', ...
       x,L(:,4),'c-', ...
       x,L(:,5),'m-');
  title('{\bf Legendre polynomails (recurrence relation)}');
  
  fig = figure('Name','Legendre polynomials (analytic values)');
  plot(x,DP_0(x),'r-', ...
       x,DP_1(x),'b-', ...
       x,DP_2(x),'g-', ...
       x,DP_3(x),'c-', ...
       x,DP_4(x),'m-');
  title('{\bf Derivatives of Legendre polynomials (analytic values)}');
  
  fig = figure('Name','Legendre polynomials (recurrence relation)');
  plot(x,DL(:,1),'r-', ...
       x,DL(:,2),'b-', ...
       x,DL(:,3),'g-', ...
       x,DL(:,4),'c-', ...
       x,DL(:,5),'m-');
  title('{\bf Derivatives of Legendre polynomails (recurrence relation)}');
  
  
  % Clear memory
  
  clear all;
  