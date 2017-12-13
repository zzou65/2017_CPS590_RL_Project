function buf = close(buf)
% CLOSE Close FIFO buffer.
%
%   CLOSE(BUF) close down a FIFO buffer BUF.
%
%   Example:
%
%   buf = open(buffer());
%   buf = close(buf);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 
  
  if(buf.open == buf.OPEN)
    buf.open = buf.CLOSED;
    buf.state = buf.READ_ONLY;
    buf.pos = 0;
    buf.size = 0;
    buf.nbytes = 0;
    try
      fclose(buf.fid);
      buf.fid = 0;
    catch
      error(ferror(fid));  
    end
    delete(buf.filename);
    buf.filename = '';
  end
  
return