function flag = isempty(buf)
% ISEMPTY Checks wheter buffer is empty.
%
%   FLAG = ISEMPTY(BUF) checks wheter the FIFO buffer BUF is empty.
%
%   Example:
%
%   isempty(buf);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  if(buf.pos < buf.nbytes)
    flag = 0;
  else
    flag = 1;
  end
      
return