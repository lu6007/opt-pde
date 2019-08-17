% --------------- section 1 ---------------
% Before use, set the parameters in section 1
% path: is the path locate the code
% init_u_tag: set initial guess got concentration, 1 is u1, 2 is u2
% init_d: initial guess for the diffusion coefficient, you can also
%         change the initial guess in later code to allow different
%         initial guess for each different coefficients
% gamma: gamma used in the algorithm
% name: decide which test problem to run
%            1-6 -- general diffusion problems with data.fun_idx = 1-6
%            'spot' -- spot diffusion
%            'layer' -- layered diffusion
%            'tensor' -- tensor diffusion
%            'tensor_cross_2' -- tensor cross diffusion
%            'mem17' -- mem17 diffusion
%            'general' -- general diffusion problem
% Example
% >> name = 1;
% >> test(name);

% Kathy 4/4/2019
% Questions:
% (1) Make >> test('general') work
% (2) **** test on H3K9 biosensor results.
% (3) split uu, vv, dd in my_newton(). 
% (4) Split the update matrix part and update diffusion coefficient part
% in cal_sub_mat(), Cal_J, and hessian_sub_mat().
% (5) Making the choice of solver an option
% (6) let v0= -A^(-1)(u0-u2)
%
% Results: 
% (1) Switched from the direct solver to the decomposition solver, it was 100x
% faster and more accurate.
% (2) Switched normalization to abs(u2-u1), converges as fast as without
% normalization.
% (3) The layered problem converges when damping was turned off in Newton's
% method
% (4) Reduced the number of input variables in to "data" in the funciton my_newton().
% (5) Changed the variable "data" from a global variable to a local
% varialbe.
% (6) Removed path and file anmes from cal_J() and hessian_sub_mat()
% (7) Made modifications and corrections in cal_sub_mat(), hessian_sub_mat(), and cal_J().
% (8) Change test_first_derivatives() to test_first_derivative().
% (9) Added a function plot_newton_step to draw the convergence history.
% (10) Removed path from cal_sub_mat()
% (11) Reduced the number of I/O variables in functions cal_sub_mat, hessian_sub_mat,
% and cal_J.
% (12) Simplified my_du() and my_newton()
% (13) Only diff_coef is scaled by a number, dt is not scaled now. 
% Updated 4/4/2019
% (14) Add an iterative solver --> GMRES solver. 
% (15) The matrices are aleady sparse
% (16) Enhance computational efficiency by using sparse matrix,
% iterative solvers and multigrid solvers. 
% (18) With damping but no normalization, sometime the algorithm does not
% converge, since the damping consant was selected assuming normalization. 
% (19) Normalize u by (u2-u1)
% (20) Improve the User interface of test.m 
% (21) Compared the speed of different solvers (direct, iterative, or preconditioned) in MATLAB. 
% 8/13/2019 and after
% (22) Switched from test_func and test_first_derivatives to
% objective functions(). 
% (23) The Mem17 case seems to be working. 
% (24) Added the opt_init_data() and load_data() functions. 
% (25) Add the codes to visualize the diffusion map. 
% (26) General unconstrained optimization works now
% (27) Checked mtool/spmatrixplot(), it seems not needed


% Copyright: Shaoying Lu and Yiwen Shi, Email: shaoying.lu@gmail.com
function test(name, varargin)
data0 = opt_init_data(name);
para.name = {'init_u_tag', 'enable_normalize', ...
    'enable_damp_newton', 'gamma', 'init_d'};
% default_value = {2, 0, 0, 1.0e-5, 10};
default_value = {data0.init_u_tag, data0.enable_normalize, ...
    data0.enable_damp_newton, data0.gamma, data0.init_d}; 
% Parse parameter and get the print_parameter() function
[init_u_tag, enable_normalize, ...
    enable_damp_newton, gamma, init_d, ~, print_parameter] = ...
    parse_parameter(para.name, default_value, varargin);
para.value = {init_u_tag, enable_normalize, ...
    enable_damp_newton, gamma, init_d};
para.format = {'%d', '%d', '%d', '%5.2e', '%6.2f'}; 
print_parameter(para);

utility_fh = utility(); % utility_function handle
update_option = data0.update_option; 
utility_fh.print_parameter('update_option', update_option);

load_data_function = data0.load_data_function; 
normalize_data_function = data0.normalize_data_function;
% ----- define the objective, jacobian and hessian functions -----
fun_idx = data0.fun_idx;
fh = objective_function; 
objective_fun = fh{fun_idx};

min_d = 0; 
max_newton_iter = data0.max_newton_iter;
max_damp_step = data0.max_damp_step;
max_outer_step = data0.max_outer_iter;

% --------------- section 2 ---------------
% This section actually runs the algorithm

% load data 
tic

data = load_data_function(name);
switch data.cell_name
    case 'general_test'
        has_constraint = 0; % No constraint
    otherwise
        has_constraint = 1; % PDE constraint
end
data.has_constraint = has_constraint; 
utility_fh.print_parameter('data.has_constraint', has_constraint);
data.gamma_d = gamma;
data = normalize_data_function(data, enable_normalize);

% ----- end -----

% # of different diffision coefficients
num_node = data.num_node;
num_para = data.num_para;

% --------------- initialization ---------------
% set the initial guess vector (u0, [v0 ,] d0)
d0 = init_d *ones(num_para, 1);
data.min_d = min_d*ones(num_para, 1);
u_len = num_node * 2 + num_para;

data.u0 = zeros(u_len, 1);
if init_u_tag == 1
  data.u0(1 : num_node) = data.u1;
elseif init_u_tag == 2
  data.u0(1 : num_node) = data.u2;
end
data.u0(num_node * 2 + 1 : end) = d0; % d0
% --------------- end ---------------

% run the actual code
data.max_newton_iter = max_newton_iter; %300;
data.max_damp_iter = enable_damp_newton * max_damp_step;
data.u_old = data.u0;
% data.u_star = cell(max_outer_iter, 1);
% data.i = cell(max_outer_iter, 1);
% data.J = cell(max_outer_iter, 1);
% data.s_hist = cell(max_outer_iter, 1);
% data.d_hist = cell(max_outer_iter, 1);
% data.u_res = cell(max_outer_iter, 1);
% data.objective = cell(max_outer_iter, 1);
% 
for outer_iter = 1 : max_outer_step
  data.outer_iter = outer_iter;
  %%%
  data = my_newton(data, objective_fun);
  %%%
  d0 = data.u_old(2*num_node+1:end);
  d_hist = data.d_hist{outer_iter};
  tmp = (d_hist(:, end)-d0)./d0;
  if has_constraint && max(abs(tmp)) < 0.01 && data.gamma_d <= 1e-5 
    break;
  end
  switch update_option
      case 1 % Update d only
          % In the outer iteration, only update gamma and d, but not u and v. 
          data.u_old(2*num_node+1:end) = data.u_new{outer_iter}(2*num_node+1:end);
      case 2 % Update u, v and d
          % Updating u and v may converge slower, but with stable convergence.
          % data.u_old(2*num_node+1:end) = data.u_new{outer_iter}(2*num_node+1:end);
          data.u_old = data.u_new{outer_iter}; 
  end
  if data.gamma_d>1e-5
      data.gamma_d = max(data.gamma_d*0.01, 1.0e-5);
      if enable_normalize
          data.gamma_d = data.gamma_d * data.xy_scale;
      end
      utility_fh.print_parameter('data.gamma_d', data.gamma_d);
  end
end

toc

% cast the cell data structure into one column vector
norm_du = data.norm_du; 
norm_dd = data.norm_dd;
du_square = (cat(1, norm_du{:})).^2;
dd_square = (cat(1, norm_dd{:})).^2;
plot_newton_step(du_square, 'Norm(du)^2');
plot_newton_step(dd_square, 'Norm(dd)^2');

u_star = data.u_new{end}; 
if size(name, 2)>1
    u_res = u_star(1:num_node);
    v_res = u_star(num_node+1 : 2*num_node);
    d_res = u_star(2*num_node+1:end); 
end

%   if size(name, 2) == 1
%       index = fun_idx; 
%       uu = data.u_res{1} + data.u2;
%       u_star = uu(end);
%       % 
%       u_min = min([uu; data.u0]);
%       u_max = max([uu; data.u0]);
%       dist = max(u_max-u_star, u_star-u_min);
%       step = 2*dist/20;
%       uuu = (u_star-dist:step: u_star+dist+step)';
%       data.x = uuu; 
%       data = f{index}(data, 0);
%       my_figure('font_size', 24, 'line_width', 3); hold on;
%       plot(data.x, data.y, 'k--', 'LineWidth', 1);
%       xlabel('u1'); ylabel('Objective');
%       
%       % 
%       data.x = [data.u0; uu];
%       data = f{index}(data, 0);
%       plot(data.x, data.y, 'b', 'LineWidth', 2);
%       plot(data.x(1), data.y(1), 'bo');
%       plot(data.x(end), data.y(end), 'r*');  
%   end

  if isfield(data, 'plot_surf') && data.plot_surf
    % d_res = d_res / data.xy_scale;  
    point = data.p'; triangle = data.tri'; edge = data.edge';
    figure; pdesurf(point', triangle', u_res);
    hold on; pdemesh(point', edge', triangle')
    colormap jet
    colorbar; view(2);
    title('u\_res');

    figure; pdesurf(point', triangle', v_res);
    hold on; pdemesh(point', edge', triangle')
    colormap jet
    colorbar; view(2);
    title('v\_res');

    figure;
    if max(size(d_res)) ~= size(triangle, 1) && max(size(d_res)) ~= data.num_node
      d_vec = zeros(1, size(triangle, 1)); 
      for i = 1 : size(triangle, 1)
        d_vec(i) = d_res(triangle(i, 4));
      end
      pdesurf(point', triangle', d_vec);
    else
      pdesurf(point', triangle', d_res');
    end
    hold on; pdemesh(point', edge', triangle')
    colormap jet; 
    colorbar; view(2);
    title('d');
  end

end
