function y = phi(x)
% implementation of a polynomial cutoff function

%   Copyright 2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


% chose from different implementations

num = 3;               


switch num
    case 1 
        y = zeros(size(x,1),1);

        case1 = x(:,2)>=abs(x(:,1));
        case2 = x(:,1)>=abs(x(:,2));
        case3 = -x(:,2)>=abs(x(:,1));
        case4 = -x(:,1)>=abs(x(:,2));

        vect = 2*x(:,2)-1;
        y(case1) = vect(case1);

        vect = 2*x(:,1)-1;
        y(case2) = vect(case2);

        vect = -2*x(:,2)-1;
        y(case3) = vect(case3);

        vect = -2*x(:,1)-1;
        y(case4) = vect(case4);
        
    case 2
        y = zeros(size(x,1),1);

        case1 = x(:,2)>=abs(x(:,1));
        case2 = x(:,1)>=abs(x(:,2));
        case3 = -x(:,2)>=abs(x(:,1));
        case4 = -x(:,1)>=abs(x(:,2));

        vect = 4*x(:,2)-3;
        y(case1) = vect(case1);

        vect = 4*x(:,1)-3;
        y(case2) = vect(case2);

        vect = -4*x(:,2)-3;
        y(case3) = vect(case3);

        vect = -4*x(:,1)-3;
        y(case4) = vect(case4);
        
        y(y<0) = 0;
        
        
    case 3
        k=7;
        
        y = zeros(size(x,1),1);

        case1 = x(:,2)>=abs(x(:,1));
        case2 = x(:,1)>=abs(x(:,2));
        case3 = -x(:,2)>=abs(x(:,1));
        case4 = -x(:,1)>=abs(x(:,2));

        vect = (2*(x(:,2)-.5)).^k;
        y(case1) = vect(case1);

        vect = (2*(x(:,1)-.5)).^k;
        y(case2) = vect(case2);

        vect = (-2*(x(:,2)+.5)).^k;
        y(case3) = vect(case3);

        vect = (-2*(x(:,1)+.5)).^k;
        y(case4) = vect(case4);
        
        
end

