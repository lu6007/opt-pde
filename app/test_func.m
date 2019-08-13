% objective function
function fh = test_func
  fh = localfunctions;
return

function data = f1(data)
  x = data.x;
  data.y = (x - 2).^2 - 5;
return

function data = f2(data)
  data.y = exp(x) - 4 .* data.x;
return

function data = f3(data)
  data.y = 2 .* (data.x - 2);
return

function data = f4(data)
  data.y = cos(x) + 1.1 * data.x;
return

function data = f5(data)
  data.y = (data.x) .^2 - 5;
return

function data = f6(data)
  data.y = (data.x) .^ (1/3);
return

% % This Jacobian, not the objective
% function [data] = f7(data)
%   x = data.x; 
% 
%   % data
%   num_node = data.num_nodes;
% 
%   u = x(1 : num_node);
%   v = x(num_node+1 : 2*num_node);
%   d = x(2*num_node+1:end);
% 
%   % assemble A and M
%   data = cal_sub_mat(data, d);
%   [A, M] = assemble_matrix(data.p, data.tri, 'diff_coef', data.diff_coef); 
%   A = A + M / data.dt;
%   data.M = M;
%   data.A = A;   
%     
%   [Ju, Jv, Jd] = cal_J(data, d, u, v);
%   y = [Ju; Jv; Jd];
%   data.y = y;
% return
