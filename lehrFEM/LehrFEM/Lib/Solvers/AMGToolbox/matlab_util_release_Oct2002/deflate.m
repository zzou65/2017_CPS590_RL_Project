%deflate		Part of block Arnoldi routine "arnoldi"
%
%usage: [q,t]=deflate(x, delta)
%
%   deflate block x such that
%   x = q*t where q is orthogonal and t has full column rank with precision
%   delta (default delta =1e-12)
%
%  Modified defaltion tolerance test.
%
%See also  arnoldi


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% These programs were prepared by the Regents of the University of
%% California at Los Alamos National Laboratory (the University) under
%% contract W-7405-ENG-36 with the U.S. Department of Energy (DOE) and
%% under Contract KC-07-01-01 with the Department of Energy(DOE),
%% Office of Science, Mathematical, Information, and Computational Sciences,
%% Applied Mathematical Sciences Program.  All rights in these programs
%% are reserved by the DOE and the University.

%% Permission is granted to the public to copy and use this software
%% without charge, provided that this Notice and the statements of 
%% authorship are reproduced on all copies. Neither the U.S. government 
%% nor the University makes any warranty, express or implied, or assumes 
%% any liability or responsibility for the use of this software.

%% AMG code written by Menno Verbeek, in collaboration with Jane Cullum
%%  and Wayne Joubert.

%% Authors Jane Cullum and Tong Zhang
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [q,t]=deflate(x, delta)

if nargin <2
  delta = 1e-12;
%  delta = 1e-7;
end

% computes permutation vector e such that q*t=x(:,e);
[q,t,e] = qr(x,0);

idebug=0;
if idebug == 1
fprintf(' in deflate permutation e equals\n');disp(e);
fprintf(' t-matrix equals\n');disp(t);
end;

for k=2:size(t,1),
% orig  if(abs(t(k,k))<delta)
  if(abs(t(k,k))/max(abs(t(1,1)),1.0) <delta)
    t=t(1:k-1,:);
    q=q(:,1:k-1);
    fprintf('deflating to  k=%5.0f\n',k-1);
    break;
  end
%fprintf(' no break in k loop in deflate k=%5.0f\n',k);
end

% ordering of new q does not matter as long as the order of t matches it

t=t(:,iorder(e));

idebug=0;
if idebug == 1
% check x - q*t is small
clear diff
qtxdiff = q*t - x;
mdiff=max(qtxdiff);
fprintf(' maximum difference between q*t -x equals %15.5e\n',mdiff);
end;




