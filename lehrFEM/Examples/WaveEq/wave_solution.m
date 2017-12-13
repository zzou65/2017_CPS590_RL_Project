function y= wave_solution(x,t,ghandle);
% exaxt solution for the wave equation
%
% Copyright 2008-2008 Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

% D'Alembert principle (Evans p. 69):
% defining the exact solution for the wave equation
% with initial conditon u given by ghandle function and v=1; 
% for reflecting boundary conditons at x = 0.0 and x=1;

a1_1 = -ghandle(-x - t);
a1_2 = ghandle(x + t);
a1_3 = -ghandle(x+t-1);

a2_1 = -ghandle( t-x);
a2_2 = ghandle(x - t);
a2_3 = -ghandle(x-t-1);

a1=a1_3;
n=length(x);
for i=1:n
    if (x(i)+t>=0 && x(i)+t<=1) a1(i)=a1_2(i); 
    else
        if (x(i)+t<0) a1(i)=a1_1(i); end
    end

    a2(i)=a2_3(i);
    if (x(i)-t>=0 && x(i)-t<=1) a2(i)=a2_2(i); 
    else
        if (x(i)-t<0) a2(i)=a2_1(i); end
    end

    y(i)=0.5*(a1(i)+a2(i));
end
y=y';