function f = indicator(x, varargin)
% implementation of the indicator function on squared and circular domains

%   Copyright 2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

domain = 1;

switch domain
    case 1
        % indicator function on the square
        
        ll = [-0. 0.5]; % left lower corner
        xspan = 0.25;
        yspan = 0.25;

        dir_1 = bitand(x(:,1) >= ll(1) , x(:,1) <=ll(1)+xspan);
        dir_2 = bitand(x(:,2) >= ll(2) , x(:,2) <=ll(2)+yspan);
        loc = bitand(dir_1,dir_2);
        
        f = zeros(size(x,1),1);
        
        f(loc) = 1;
        
    case 2
        % indicator function on a circle

        m = [0 0.5];
        r = 0.25;

        m = [m(1)*ones(size(x,1),1) m(2)*ones(size(x,1),1)];

        rsquared = sum((x - m).^2,2);

        loc = rsquared <= r^2;

        f = zeros(size(x,1),1);
        
end

f(loc) = 1;