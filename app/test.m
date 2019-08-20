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

% Kathy 8/19/2019
% Questions:
% (1) Make >> test('general') work
% (2) **** test on H3K9 biosensor results.
% (3) Making the choice of solver an option
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
% (28) split uu, vv, dd in my_newton(). 

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
u_len = num_node;

data.u0 = zeros(u_len, 1);
if init_u_tag == 1
  data.u0 = data.u1;
elseif init_u_tag == 2
  data.u0 = data.u2;
end
data.v0 = zeros(u_len, 1); 
data.d0 = d0; % d0
% --------------- end ---------------

% run the actual code
data.max_newton_iter = max_newton_iter; %300;
data.max_damp_iter = enable_damp_newton * max_damp_step;
data.u_old = data.u0; data.v_old = data.v0; data.d_old = data.d0; 
for outer_iter = 1 : max_outer_step
  data.outer_iter = outer_iter;
  %%%
  data = my_newton(data, objective_fun);
  %%%
  if ~has_constraint
      continue;
  end
  
  % has constraint
  tmp = (data.d_new-data.d_old)./data.d_old;
  data.d_old = data.d_new; 
  if update_option == 2 % Update u and v
      % Updating u and v may converge slower, but with stable convergence.
      data.u_old = data.u_new;
      data.v_old = data.v_new; 
  end
  if max(abs(tmp)) < 0.01 && data.gamma_d <= 1e-5 
      break;
  end
  
  %
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

if isfield(data, 'draw_surface') && data.draw_surface
    dfun = diffusion_function(); % require diffusion_analysis
    surf.p = data.p; surf.tri = data.tri; surf.u = data.u_new';
    surf.node = data.p; surf.edge = data.edge;
    figure; 
    dfun.draw_surface_with_mesh(surf, surf); 
    caxis([0 6.5]); title('u\_star'); 

    surf.u = data.v_new';
    figure; dfun.draw_surface_with_mesh(surf, surf); 
    title('v\_star');

    d_new = data.d_new * data.xy_scale; 
    triangle = data.tri'; 
    if size(d_new, 1) ~= size(triangle, 1) && size(d_new, 1) ~= data.num_node
      d_vec = zeros(1, size(triangle, 1)); 
      for i = 1 : size(triangle, 1)
        d_vec(i) = d_new(triangle(i, 4));
      end
      surf.u = d_vec';
    else
      surf.u = d_new;
    end
    figure; dfun.draw_surface_with_mesh(surf, surf); 
    caxis([0 250]); title('d\_star');
end

return
