%ainv		Compute approximate factorization of inverse of SPD matrix
%
%usage : [Z,D]=ainv(A,imind,dtol)
%
% Compute approximate factorization of inverse of
%  a positive definite symmetric matrix A
%
% Procedure constructs A-orthogonal set of vectors Z
%  from the starting vectors Z0 which by default
%  is set to the Identity matrix
%
% Input:
%   A = positive definite symmetric matrix
%   imind = 1; apply minimum degree reordering to A; 
%   dtol = drop tolerance applied to z-vectors
%
% Output:
%   Z = upper triangular matrix: columns are z-vectors generated
%        after drop tolerance is applied


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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [Z,D] = ainv(A,imind,dtol)



% Check matrix is square
   [nA,nAc]= size(A);
   if (nA ~= nAc) 
     disp(' A is not Square: Terminate');
     return;
%
% Proceed with constructing the preconditioner
   else
%
% Default Z starting block: Identity matrix
   if (nargin <= 4);  Z0 = speye(size(A)); end;

% Requires symmetric A

%disp(' Flag imind');disp(imind);
%disp(' Drop Tolerance dtol');disp(dtol);

   
% Benzi typically applies reordering prior to 
%  constructing a preconditioner

    if imind ==  1

% Reorder using Minimum degree 
% Symmmd = symmetric minimum degree permutation

     clear P;
     P = symamd(A);

% Apply permutation to A

     A=A(P,P);
    end;

% Construct Preconditioner
% Left looking:

% Use z_i'*A*z_k: 

% Generate Z and D such that Z^T * A * Z = D
%  where D is diagonal, Z is block upper triangular

% In exact arithmetic and no dropping, A^{-1} = Z * D^{-1} * Z^T

% Loop on column in Z
% Z and D are matlab sparse matrices
% Z0 contains starting Z matrix in sparse matlab format.


% Temp holds A*z_j for reuse

%    D = sparse(nA,nA);
%    Z = sparse(nA,nA);
%    Temp = sparse(nA,nA);

% initialize number of z vector
    i = 1;
% initialize number of blocks with more than 1 vector
    nb = 0;
% initialize number of vectors in a given block
    isb = 1;
% Generate the z-vectors
%
  while i <= nA
%
% this assumes Z0 specified
% we could consider other choices for Z0 than the identity matrix.

     Z(1:nA,i)= Z0(:,i);
%
%   A-orthogonalize new Z vector w.r.t existing Z-vectors
%

   if i ~= 1
   if 1
%'new'
     % Do 2 steps of Ordinary GS
     J=1:i-1;
     alpha=Temp(:,J)'*Z(:,i);
     alpha=alpha./diag(D(J,J));
     Z(:,i)=Z(:,i)-Z(:,J)*alpha;
     alpha=Temp(:,J)'*Z(:,i);
     alpha=alpha./diag(D(J,J));
     Z(:,i)=Z(:,i)-Z(:,J)*alpha;

   else

     % Do MGS
       Zi=full(Z(:,i));
       for j=1:i-1
         % Use A*z_j :

         %cij = Z(:,j)'*A*Z(:,i);
         %cij = Temp(:,j)'*Z(:,i);
         cij = Zi'*Temp(:,j);

         %Z(:,i) = Z(:,i) - (cij/D(j,j))*Z(:,j);
         Zi = Zi - (cij/D(j,j))*Z(:,j);
       end;
       Z(:,i)=Zi;

   end
   end;

%

% Sparsify  new vector
   is  = i - isb;
%
   for k=1:isb;
     is = is + 1;  
     inz=find(Z(:,is));
     nz =length(inz);
     for m=1:nz
% if entry small and not a diagonal entry
      if (abs(Z(inz(m),is)) < dtol ) & (inz(m) ~= is)
% drop term
        Z(inz(m),is) = 0;
      end;
     end;
   end;

% Save A*Z(:,i)
   Temp(1:nA,i)=A*Z(:,i);

%   Compute new diagonal term
%
%     D(i,i)= Z(:,i)'*A*Z(:,i);
      D(i,i)= Temp(:,i)'*Z(:,i);


% increment vector number
   i = i + 1;
% reset number of vectors in block
   isb = 1;
% end of while
 end;

   if (imind == 1)

%  remove the symmmd permutation
%  this assume A_permuted = P*A*Pinv
%  ie replace Z by Z(pinv,I) i.e. 
%   inv(A(p,p)) = P*inv(A)*Pinv = Z*inv(D)*Z^T
%   implies 
%       inv(A)= Pinv*Z *inv(D)*Z^T*P

% compute Pinv
    [dummy, Pinv]=sort(P);
    Z= Z(Pinv,:);
  
   end;

% end of if matrix is not square
  end;





