%% NCC_TRIANGLE_TEST is the main program for the NCC_TRIANGLE sample code.
%
%  Modified:
%
%    12 December 2006
%
%  Author:
%
%    John Burkardt
%
  fprintf ( 1, '\n' );
  timestamp;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'NCC_TRIANGLE_PRB:\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Demonstrate routines in the NCC_TRIANGLE library.\n' );

  ncc_triangle_test01;
  ncc_triangle_test02;
  ncc_triangle_test03;
  ncc_triangle_test04;
  ncc_triangle_test05;
  ncc_triangle_test06;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'NCC_TRIANGLE_PRB:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );

  fprintf ( 1, '\n' );
  timestamp;
