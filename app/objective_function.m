% function fh = objective_function
%   fh = localfunctions;
% return
% Objective functions,
% func_flab = 0 --- The objective function
% func_flag = 1 --- First derivative, or Jacobian
% func_flag = 2 --- Second derivative, or Hessian

% Copyright: Shaoying Lu and Yiwen Shi, Email: shaoying.lu@gmail.com
function fh = objective_function
  fh = localfunctions;
return

function data = f1(data, func_flag)
  x = data.x;
  if func_flag == 0
    % objective function
    data.y = (x - 2)^2 + 5;
  elseif func_flag == 1 
    % jacobian function
    data.y = 2 * (x - 2);
  elseif func_flag == 2
    % hessian function
    data.y = 2;
  end

return

function data = f2(data, func_flag)
  x = data.x;
  if func_flag == 0
    % data.y = exp(x) - 2 * x.^2;
    data.y = exp(x)-2*x^2;
  elseif func_flag == 1
    data.y = exp(x) - 4 * x;
  elseif func_flag == 2
    % data.y = diag(exp(x) - 4);
    data.y = exp(x)-4; 
  end
return

function data = f3(data, func_flag)
  x = data.x;
  if func_flag == 0
    data.y = sum((x - 2).^2);
  elseif func_flag == 1
    data.y = 2 * (x - 2);
  elseif func_flag == 2
    data.y = diag(2 * ones(length(x), 1));
  end
return

function data = f4(data, func_flag)
  x = data.x;
  if func_flag == 0
    data.y = sum(sin(x) + 1.1/2 * x.^2);
  elseif func_flag == 1
    data.y = cos(x) + 1.1 * x;
  elseif func_flag == 2
    data.y = -sin(x) + 1.1;
  end
return

function data = f5(data, func_flag)
  x = data.x;
  if func_flag == 0
    data.y = sum(x.^3 / 3 - 5 * x);
  elseif func_flag == 1
    data.y = (x) .^2 - 5*ones(length(x), 1);
  elseif func_flag == 2
    data.y = 2 * diag(x);
  end
return

% f6 is an example of non-convergence
function data = f6(data, func_flag) 
  x = data.x;
  if func_flag == 0
    % data.y = 3/4 * nthroot(x,3).^4; % sum of 3/4 * x.^(4/3)
    data.y = sum(3/4 * x.^(4/3)); 
  elseif func_flag == 1
    % Use the nthroot function to return a real number
    data.y = nthroot(x, 3); % x.^(1/3)
    % data.y = x.^(1/3);
  elseif func_flag == 2
    data.y = 1/3 * diag(nthroot(x, 3) .^ (-2)); 
    % data.y = 1/3 * diag(x.(-2/3)); 
  end
return

function [data] = f7(data, func_flag)
  x = data.x;
  u = x.u; v = x.v; d = x.d; 
  
  if func_flag == 0 % objective
    gamma_d = data.gamma_d; 
    norm_u = (u - data.u2)' * data.M * (u - data.u2); 
    AM = data.A; % data.A = A+M/data.dt
    constraint = u' * AM * v + data.u1' * data.M * v / data.dt; 
    norm_d = 0;
    if min(abs(data.d_old))<=eps
        fprintf('\n Function objective_function(): some abs(d0)<=eps, do not normalize. \n');
    end
    for iter = 1 : data.num_para
      tmp_M_sub = data.M_sub{iter};
      dd0 = data.d_old(iter);
      ddd = d(iter)-dd0;
      if abs(dd0) > eps
        norm_d = norm_d + gamma_d * (ddd^2) * sum(tmp_M_sub(:)) /(dd0^2);
      else
        norm_d = norm_d + gamma_d * (ddd^2) * sum(tmp_M_sub(:));
      end
      clear tmp_M_sub; 
    end
    % 
    y = norm_u + norm_d + constraint;
    
  elseif func_flag == 1 % update matrices and calculate first derivative
    data = cal_sub_mat(data, d);
    % Assemble and update the matrics and submatrices
    % assemble A and M
    data.assemble = 1; 
    if ~isfield(data, 'assemble') || data.assemble
      [A, M] = assemble_matrix(data.p, data.tri, 'diff_coef', data.diff_coef); 
      A = A + M / data.dt;
      data.M = M;
      data.A = A;   % changing matrix A here, very confusing. when was data.A first assembled? 
    end
    %
    [Ju, Jv, Jd] = cal_J(data, d, u, v);
    y = [Ju; Jv; Jd];
  elseif func_flag == 2 % calculate second derivative
    % Call objective_function with fun_flag = 1 first to assemble A_sub and M_sub
    % Then update Cu, Cv, and G
    data = hessian_sub_mat(data, u, v);
    Cu = data.Cu;
    Cv = data.Cv;
    G = data.G;
    if max(G(:)) == Inf
      fprintf('\nWarning --- objective_function(): max(G(:)) == Inf\n\n');
      pause
    end
    A = data.A;
    M = data.M;
    O = zeros(size(A));
    H = [M ,  A,   Cu;
         A,   O,   Cv;
         Cu', Cv', G ];

    y = H;
    clear Cu Cv G H; 
  end
  data.y = y;
return

% Not sure about this purpose, need to be tested. 8/12/2019 Kathy
function data = f8(data, func_flag)
    data = f7(data, func_flag);
    num_node = data.num_node;
    x = data.x;
    d = x(2*num_node+1:end);
    gamma2 = 1; 
    % Calculate the H1 half norm of d and its derivatives. 
    switch func_flag
        case 0 % objective
            dh1_half = 0; 
    end
    data.y = data.y + gamma2*dh1_half;
return

% function data = objective_diffusion_h1norm(data, func_flag)
%   x = data.x;
%   num_node = data.node;
% 
%   u = x(1 : num_node);
%   v = x(num_node+1 : 2*num_node);
%   d = x(2*num_node+1:end);
%   
%   p = data.p; 
%   tri = data.tri;
% 
%   if func_flag == 0
% 
%     objective_tmp = (u - data.u2)' * data.M * (u - data.u2) ...
%                   + u' * data.A * v - data.u1' * data.M * v / data.dt;
%                         % + data.u1' * data.M * u_tmp(num_node+1:2*num_node) / data.dt;
% 
%     tmp_obj = 0;
%     if strcmp(data.cell_name, 'layered_diffusion_general_h1')
%       delta_d = (d-data.d_old) ./ data.d_old;
%       [A_tmp, M_tmp] = assemble_matrix(p', tri');
% %       tmp_obj = data.gamma_d * delta_d' * data.M * delta_d...
% %          + data.gamma_d * delta_d' * (data.A - data.M / data.dt) * delta_d;
%       tmp_obj = data.gamma_d * delta_d' * M_tmp * delta_d...
%         + data.gamma_d * delta_d' * A_tmp * delta_d;
%     else
%       for iter = 1 : data.num_para
%         tmp_M_sub = data.M_sub{iter};
%         d0 = data.d_old(iter);
%         delta_d = d(iter)-d0;
%         tmp_obj = tmp_obj + data.gamma_d * (delta_d^2) * sum(tmp_M_sub(:)) /(d0^2);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tmp_A_sub = data.A_sub{iter};
%         % tmp_obj = tmp_obj - u_tmp(2*num_node+iter)^2 * sum(sum(data.A_sub{iter}));
%         % tmp_obj = tmp_obj + u_tmp(2*num_node+iter)^2 * sum(sum(data.A_sub{iter}));
%         % tmp_obj = tmp_obj + u_tmp(2*num_node+iter)^2 * sum(sum(data.A_sub{iter})) / data.gamma_d;
%         % tmp_obj = tmp_obj + u_tmp(2*num_node+iter)^2 * sum(sum(abs(data.A_sub{iter})));
%         % tmp_obj = tmp_obj + u_tmp(2*num_node+iter)^3 * sum(sum(data.A_sub{iter}));
%         % tmp_obj = tmp_obj + u_tmp(2*num_node+iter)^3 * sum(sum(abs(data.A_sub{iter})));
%         tmp_obj = tmp_obj + data.gamma_d * (delta_d^2) * sum(tmp_A_sub(:)) / (d0^2);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       end
%     end
%     objective_tmp = objective_tmp + tmp_obj;
%     data.y = objective_tmp;
%   elseif func_flag == 1
%     % data
%     % assemble A and M 
%     if ~isfield(data, 'assemble') || data.assemble
%       data = cal_sub_mat(data, d);
% 
%       [A, M] = assemble_matrix(p', tri', 'diff_coef', data.diff_coef');
%       A = A + M / data.dt;
%       data.M = M;
%       data.A = A;
%     end
% 
%     [Ju, Jv, Jd] = cal_J(data, d, u, v);
%     y = [Ju; Jv; Jd];
%     data.y = y;
%   elseif func_flag == 2
%     % data = hessian_sub_mat(data, d, u, v);
%     data = hessian_sub_mat(data, u, v);
% 
%     Cu = data.Cu;
%     Cv = data.Cv;
%     G = data.G;
%     if max(G(:)) == Inf
%       pause
%     end
%     A = data.A;
%     M = data.M;
%     fillin = zeros(size(A));
%     H = [M ,   A,      Cu;
%          A,   fillin, Cv;
%          Cu', Cv',    G ];
% 
%     data.y = H;
%   end
% 
% return