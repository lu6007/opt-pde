function fh = test_first_derivative
  fh = localfunctions;
return

function data = d1(data)
  x = data.x;
  data.y = 2 .* (x - 2);
return

function data = d2(data)
  data.y = diag(exp(data.x) - 4);
return

function data = d3(data)
  data.y = diag(2 * ones(length(data.x), 1));
return

function data = d4(data)
  data.y = sin(data.x) + 1.1;
return

function data = d5(data)
  data.y = 2 .* (data.x);
return

function data = d6(data)
  data.y = 1/3 .* (data.x) .^ (-2/3);
return

function data = d7(data) % This is Hessian, not the first derivative. 
  x = data.x;
  num_node = data.num_nodes; 

  u = x(1 : num_node);
  v = x(num_node + 1 : 2 * num_node);
  % d = x(2 * num_node + 1 : end);

%   Call test_func() first to assemble A_sub, M_sub, A and M
% Update Cu, Cv, and G
  data = hessian_sub_mat(data, u, v);

  Cu = data.Cu;
  Cv = data.Cv; 
  G = data.G; 
  A = data.A; M = data.M;
  fillin = zeros(size(A));
  H = [M,   A,      Cu;
       A,   fillin, Cv;
       Cu', Cv',    G ];

  data.y = H;
  clear Cu Cv G H;
return
