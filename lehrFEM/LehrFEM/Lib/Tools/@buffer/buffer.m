function buf = buffer(arg)
% BUFFER Open up a new FIFO buffer.
%
%   BUF = BUFFER() opens up a new FIFO buffer.
%
%   Example:
%
%   buf = buffer();

%   Copyright 2006-2006 Patrick Meuryfclose(buf.fid);
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

  % Initialize static buffer constants
  
  buf.OPEN = 0;        % Buffer is open
  buf.CLOSED = 1;      % Buffer is closed
  buf.READ_ONLY = 2;   % Buffer is read only
  buf.WRITE_ONLY = 3;  % Buffer is write only
  buf.REAL = 4;        % Real valued data
  buf.COMP = 5;        % Complex valued data
  
  % Initialize private buffer variables
  
  buf.open = buf.CLOSED;      % Buffer file state indicator
  buf.state = buf.READ_ONLY;  % Buffer state indicator
  buf.filename = '';          % Filename of buffer file
  buf.fid = 0;                % File identifier
  buf.type = 0;               % Data type flag
  buf.pos = 0;                % Position indicator [bytes from BOF]
  buf.size = 0;               % Number of elements
  buf.nbytes = 0;             % Size of buffer file [bytes from BOF]
     
  buf = class(buf,'buffer');
  
return