% --------------- section 1 ---------------
% Before use, set the parameters in section 1
% path: is the path locate the code
% init_u_tag: set initial guess got concentration, 1 is u1, 2 is u2
% init_d: initial guess for the diffusion coefficient, you can also
%         change the initial guess in later code to allow different
%         initial guess for each different coefficients
% gamma: gamma used in the algorithm
% test_case: decide which test case to run
%            1 -- spot diffusion
%            2 -- layered diffusion
%            3 -- tensor diffusion
%            4 -- tensor cross  diffusion
%            5 -- mem17 diffusion

% Kathy 4/4/2019
% Questions:
% (1) How to run a general test function for a simple objective funtion to
% demonstrate the simple normalization and damping Newton's method works?
% --> how to initialize for a general function? (Yiwen)
% (3) Split the update matrix part and update diffusion coefficient part
% in cal_sub_mat(), Cal_J, and hessian_sub_mat().
% (5) **** test on H3K9 biosensor results.
% (6) mtool/spmatrixplot()
% (7) Making the choice of solver an option
% (8) let v0= -A^(-1)(u0-u2)
% (9) Test case 7: not convergent yet. 

%
% Results
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
% 

% Copyright: Shaoying Lu and Yiwen Shi, Email: shaoying.lu@gmail.com
function test(varargin)
para.name = {'test_case', 'init_u_tag', 'enable_normalize', ...
    'enable_damp_newton', 'gamma', 'init_d'};
default_value = {1, 2, 0, 0, 1.0e-5, 10};
% Parse parameter and get the print_parameter() function
[test_case, init_u_tag, enable_normalize, ...
    enable_damp_newton, gamma, init_d, ~, print_parameter] = ...
    parse_parameter(para.name, default_value, varargin);
para.value = {test_case, init_u_tag, enable_normalize, ...
    enable_damp_newton, gamma, init_d};
para.format = {'%d', '%d', '%d', '%d', '%5.2e', '%6.2f'}; 
print_parameter(para);

utility_fh = utility(); % utility_function handle
update_option = 2; % In the outer iterations, update d only, works better for accurate data
% for mem17 use 2
% update_option = 2; % Update all three variables, u, v and d, works better
% for noisy data. 
utility_fh.print_parameter('update_option', update_option);

% test_case = 1; % 1 -- spot diffusion
% init_u_tag = 1; % 1 -- u0=u1; 2 -- u0 = u2
% init_d = 30; % both 10 and 30 converges
min_d = 0; 
max_newton_iter = 10;
max_damp_step = 10;
max_outer_step = 10;

% --------------- section 2 ---------------
% This section actually runs the algorithm

% load data and make necessary modifications
tic
% data
pa = '../data/';
switch test_case
  case 1
    data_file = strcat(pa, 'spot_diffusion.mat');
    data = load(data_file);
    data.tri(4, data.tri(4, :) == 3) = 2;
    data.tri(4, data.tri(4, :) == 4) = 2;
    data.tri(4, data.tri(4, :) == 5) = 2;
    data.cell_name = 'spot_diffusion';
  case 2
    data_file = strcat(pa, 'layered_diffusion.mat');
    data = load(data_file);
    data.cell_name = 'layered_diffusion';
  case 3
    data_file = strcat(pa, 'tensor_diffusion.mat');
    data = load(data_file);
    data.cell_name = 'tensor_diffusion';
    % data.KK = 0; 
  case 4
    data_file = strcat(pa, 'tensor_cross_2.mat');
    data = load(data_file);
    data.tri(4, data.tri(4, :) == 3) = 2;
    data.tri(4, data.tri(4, :) == 4) = 3;
    data.tri(4, data.tri(4, :) == 5) = 3;
    data.tri(4, data.tri(4, :) == 6) = 4;
    data.cell_name = 'tensor_cross_2';
    % data.num_para = 3;
  case 5
    data_file = strcat(pa, 'mem17.mat');
    data = load(data_file);
    data.cell_name = 'mem17';
    data.u1 = data.u(:, 1);
    data.u2 = data.u(:, 2);
    data.p = data.p_image;
    data.dt = data.dt * 60; % Convert min to sec
    mag = 98;
    data.p = scale_by_magnification(data.p, mag); % requires fluocell
    utility_fh.print_parameter('magnification', mag);
  case 6
    data_file = strcat(pa, 'WTH3K9_4.mat');
    data = load(data_file);
    data.cell_name = 'WTH3K9';
    data.u1 = data.u(:, 4);
    data.u2 = data.u(:, 5);
    data.p = data.p_image;
    data.dt = data.dt * 60;
    mag = 99; 
    data.p = scale_by_magnification(data.p, mag); % requires fluocell
    utility_fh.print_parameter('magnification', mag);
    % Use the computer simulation program to mark regions 1 and 2
    data.tri(4, data.tri(4, :) == 1) = 2;
    data.tri(4, data.tri(4, :) == 3) = 2;
    data.tri(4, data.tri(4, :) == 4) = 2;
    data.tri(4, data.tri(4, :) == 5) = 1;
  case 7 % Testing case not working yet. 8/7/2019
    data_file = strcat(pa, 'layered_diffusion_general_5_refined.mat');
    data = load(data_file);
    for i = 1 : size(data.tri, 2)
      data.tri(4, i) = i;
    end
    data.cell_name = 'layered_diffusion_general';
    data.p = data.p_image; 
end
data.path = pa;
data.num_para = size(unique(data.tri(4, :)), 2);
if test_case == 4
    data.num_para = 3;
end

data.gamma_d = gamma;

if enable_normalize
    % --------------- normalization ---------------
    % zmax = max(max(data.u1),max(data.u2));
    % zmin = min(min(data.u1),min(data.u2));
    abs_u_diff = abs(data.u2-data.u1);
    umax = max(abs_u_diff);
    umin = min(abs_u_diff);
    xmax = max(data.p(1, :));
    xmin = min(data.p(1, :));
    ymax = max(data.p(2, :));
    ymin = min(data.p(2, :));
    xmid  = (xmax+xmin)/2.0;
    ymid  = (ymax+ymin)/2.0;
    xy_scale = 1.0 / max(xmax-xmin,ymax-ymin);
    u_scale = 1.0 / (umax-umin);

    p1 = xy_scale * (data.p(1, :) - xmid);
    p2 = xy_scale * (data.p(2, :) - ymid);
%     data.u1 = (data.u1 - umin)*u_scale;
%     data.u2 = (data.u2 - umin)*u_scale;
    data.u1 = data.u1 * u_scale;
    data.u2 = data.u2 * u_scale; 
    data.p = [p1;p2];
    data.dt = data.dt * (xy_scale^2);
    data.gamma_d = data.gamma_d * xy_scale;
    % The diffusion coefficients do not need to be scaled. 
    
    % --------------- end ---------------
else
    xy_scale = 1;
end
data.xy_scale = xy_scale;

% ----- define the jacobian function and hessian -----
% fd = test_first_derivative; % Hessien
% f = test_func; % Jacobian
fh = objective_function; 
fun_idx = 7;
objective_fun = fh{fun_idx};
% ----- end -----

% # of different diffision coefficients
num_node = data.num_nodes;
num_para = data.num_para;

% --------------- initialization ---------------
% set the initial guess vector (u0, v0, d0)
u_len = num_node * 2 + num_para;
d0 = init_d * ones(num_para, 1);
data.min_d = min_d*ones(num_para, 1); 

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
% 
for outer_iter = 1 : max_outer_step
  data.outer_iter = outer_iter;
  data = my_newton(data, objective_fun);
  d0 = data.u_old(2*num_node+1:end);
  d_hist = data.d_hist{outer_iter};
  tmp = (d_hist(:, end)-d0)./d0;
  if max(abs(tmp)) < 0.01 && data.gamma_d <= 1e-5 
    break;
  end
  switch update_option
      case 1 %Update d only
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
          data.gamma_d = data.gamma_d*xy_scale;
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

% d_res = d_res / xy_scale 
end
