function n = size(buf)
% SIZE Returns the current size of a FIFO buffer.
%
%   N = SIZE(BUF) returns the current number of elements stored inside the
%   FIFO buffer BUF.
%
%   Example:
%
%   n = size(buf);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  n = buf.size;

return