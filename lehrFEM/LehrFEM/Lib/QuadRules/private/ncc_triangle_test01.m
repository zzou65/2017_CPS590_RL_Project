%% TEST01 tests NCC_TRIANGLE_RULE_NUM, NCC_TRIANGLE_DEGREE, NCC_TRIANGLE_ORDER_NUM.
%
%  Modified:
%
%    29 January 2007
%
%  Author:
%
%    John Burkardt
%
  clear

  fprintf ( 1, '\n' );
  fprintf ( 1, 'TEST01\n' );
  fprintf ( 1, '  NCC_TRIANGLE_RULE_NUM returns the number of rules;\n' );
  fprintf ( 1, '  NCC_TRIANGLE_DEGREE returns the degree of a rule;\n' );
  fprintf ( 1, '  NCC_TRIANGLE_ORDER_NUM returns the order of a rule.\n' );

  rule_num = ncc_triangle_rule_num ( 'DUMMY' );

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Number of available rules = %d\n', rule_num );
  fprintf ( 1, '\n' );
  fprintf ( 1, '      Rule    Degree     Order\n' );
  fprintf ( 1, '\n' );

  for rule = 1 : rule_num
    order_num = ncc_triangle_order_num ( rule );
    degree = ncc_triangle_degree ( rule );
    fprintf ( 1, '  %8d  %8d  %8d\n', rule, degree, order_num );
  end
