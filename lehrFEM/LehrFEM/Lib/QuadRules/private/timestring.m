function s = timestring

%% TIMESTRING returns a string containing the current YMDHMS date.
%
%  Modified:
%
%    07 August 2003
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Output, string S, a string containing the current YMDHMS date.
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );

  
