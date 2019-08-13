% u0: initial guess, including u0, v0, and d0
% s: damping step
% jacobian_fun: jacobian of objective function
% hessian_fun: jocabian matrix of jacobian function, aka hessian of objective function
% max_newton_iter: maximum number of newton iterations
% max_damp_iter: maximum number of damping steps
% data should provide meshing information, and initial guess

function data = my_newton(data, objective_fun)

%   % function handles
%   fh = objective_fun; 

  % data
  max_newton_iter = data.max_newton_iter;
  max_damp_iter = data.max_damp_iter;
  num_node = data.num_nodes;
  num_para = data.num_para; 
  gamma_d = data.gamma_d; 
  min_d = data.min_d; 
  u_old = data.u_old;
  tol = 1.0e2 * eps;

  objective = zeros(max_newton_iter+1, 1);
  J  = zeros(max_newton_iter+1, 1);
  s_hist = zeros(max_newton_iter+1, 1);
  norm_du = zeros(max_newton_iter+1, 1);
  d_hist = zeros(num_para, max_newton_iter+1);
  norm_dd = zeros(max_newton_iter+1, 1);
  % data.Cv_bar and data.d0 are only used inside the my_newton function, 
  % but not outside. 
  data.Cv_bar = zeros(num_node, num_para);
  data.d0 = data.u_old(2 * num_node+1: end); 

  % Newton iterations
  % Objective Function, Norm of Jacobian, Damping Parameter
  % Residual, Diffusion Coefficients
  fprintf('Obj Func \t Norm J \t Damping \t ');
  fprintf('Residual \t Diffusion Coefficients\n');
  %
  norm_du(1) = norm(u_old(1:num_node) - data.u2);
  norm_dd(1) = 0; 
  d_hist(:, 1) = data.d0;
  for i = 1: max_newton_iter
    data.x = u_old;
    data = objective_fun(data, 1);
    jacobian_old = data.y;
    if i == 1
        J(i) = norm(jacobian_old);
        data = objective_fun(data, 0);
        objective(i) = data.y;
    end
    % data = hessian_fun(data);
    data = objective_fun(data, 2); 
    hessian_old = data.y;

    if J(i) < tol
        u_new = u_old; 
        break;
    end
    
%     %%%%%%%%%% Sovlvers %%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Direct solver
%     du = - hessian_old \ jacobian_old;
%     %
%     % GMRES with ILU preconditioner
%     [L, U] = ilu(hessian_ old, struct('type', 'ilutp', 'droptol', 1e-5));
%     du = - gmres(hessian_old, jacobian_old, [], 1.0e-6, 2000, L, U);
%     %
    % OPT-PDE solver 
    data = my_du(data, jacobian_old);
    du = [data.delta_u; data.delta_v; data.delta_d]; 
    %
    
    if ~max_damp_iter % no damping
        ss = 1; 
        u_new = u_old + ss * du;
        data.x = u_new;
        data = objective_fun(data, 1);
        jacobian_new = data.y; 
    else
        % Damping
        alpha0 = norm(hessian_old * du + jacobian_old) / J(i);
        delta = (1 - alpha0) / 2; % set delta: 0 < delta < 1 - alpha0
        % KK = data.KK;
        KK = 0;
        for damp_iter = 1:max_damp_iter
            ss = 1.0/(1+KK * J(i));
            %%%%%%%%%%%%%%%%%%%%%
            u_new = u_old + ss * du;
            data.x = u_new;
            % data = jacobian_fun(data);
            data = objective_fun(data, 1);
            jacobian_new = data.y;
            %%%%%%%%%%%%%%%%%%%%%
            if (1 - norm(jacobian_new)/J(i))/ss < delta
              if KK == 0
                KK = 1;
              else
                KK = 10 * KK;
              end
            else
                data.KK = KK / 10;
                break;
            end
        end % for damp_iter = 1:max_damp_iter
        % fprintf('# of damping steps at step %d: %d \n', i, damp_iter);
    end % if ~max_damp_iter % no damping
    
    % constraint dd >= min_d
    % This step is necessary in some cases. 
    % This bound is only optionally triggered. Need to print something
    % here. 
    u_new(2*num_node+1: end) = max(u_new(2*num_node+1: end),min_d);
    
    % calculate objective function and Jacobian
    uu = u_new(1:num_node);
    % vv = u_new(num_node+1:2*num_node); 
    dd = u_new(2*num_node+1: end); 
    data.x = u_new;
    data = objective_fun(data, 0);
    objective(i+1) = data.y;
    J(i+1) = norm(jacobian_new);

    d_hist(:, i+1) = dd;
    s_hist(i+1) = ss;
    norm_du(i+1) = norm(uu - data.u2);
    norm_dd(i+1) = norm(dd - d_hist(:,i))/norm(data.d0);

  % From left to right: Objective Function, Norm of Jacobian, Damping Parameter
  % Residual, Diffusion Coefficients
    fprintf('%8.3f \t %10.5f \t %5.2e \t %5.2e \t %s \n', ...
        objective(i), J(i), s_hist(i), norm_du(i), num2str(u_new(2*num_node+1:end)'));

    % convergence test
    tmp = u_old(2*num_node+1: end) - dd;
    tmp = tmp ./ u_old(2*num_node+1: end);

    if J(i+1) < 1e-06 && max(abs(tmp)) < 0.01
      break;
    else
      clear u_old jacobian_old;
      u_old = u_new;
      data.x = u_old;
      % data = hessian_fun(data);
      data = objective_fun(data, 2); % Hessian
    end % if norm(jacobian_new) < 1e-06 && max(abs(tmp)) < 0.01
  end % for i = 1 : max_newton_iter

  outer_iter = data.outer_iter;
  num_newton_iter = i+1;
  data.u_new{outer_iter} = u_new;
  data.i{outer_iter} = i;
  data.J{outer_iter} = J;
  data.s_hist{outer_iter} = s_hist(1:num_newton_iter);
  data.d_hist{outer_iter} = d_hist(:, 1:num_newton_iter);
  data.norm_du{outer_iter} = norm_du(1:num_newton_iter);
  data.norm_dd{outer_iter} = norm_dd(1:num_newton_iter)*sqrt(num_node*gamma_d);
  data.objective{outer_iter} = objective(1:num_newton_iter);

  fprintf('Number of newton steps: %d \n', i);
  fprintf('######################## \n\n');
  
return

