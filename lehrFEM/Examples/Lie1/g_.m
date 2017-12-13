function g= g_(x,varargin)
% F_LSHAP Right hand side source.
%
%   Boundary values
%   Example:
%
%   g = g_([0 0]);
%
%   Copyright 2005-2006 Patrick Meury & Kah Ling Sia & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Compute boundary term

  s=size(x,1);
  g = zeros(size(x,1),2);
  
  for i = 1:s
      if ( x(i,1) == (-1))
          g(i,:)=[1 1];
      end    
      if ( x(i,2) == (1))
         g(i,:)=[0 0];
      end
      if ( x(i,1) == (1))
       g(i,:)=[0 0];   
      end
      if ( x(i,2) == (-1))
       if (x(i,1) <=1/3)
           g(i,:)=[1 1];
       else g(i,:)=[0 0];
       end    
      end
  end
% 
% 
%
%  for i = 1:s
%       if ( x(i,1) == (0) || x(i,2) == 1)
%           g(i)=g(i)+1;
%       end
%   end 
%   
%   for i = 1:s
%       if ( x(i,1) == (-1))
%           g(i)=g(i)-1+x(i,2);
%       else    
%         if ( x(i,2) == (1))
%           g(i)=g(i)+1+x(i,1);
%         end
%       end
%  end

%     if ( x(i,1) == (-1) )
%         g(i)=g(i)+1;
%     end
%  end
return