% u0: initial guess, including u0, v0, and d0
% s: damping step
% jacobian_fun: jacobian of objective function
% hessian_fun: jocabian matrix of jacobian function, aka hessian of objective function
% max_newton_iter: maximum number of newton iterations
% max_damp_iter: maximum number of damping steps
% data should provide meshing information, and initial guess

% Copyright: Shaoying Lu and Yiwen Shi, Email: shaoying.lu@gmail.com
function data = my_newton(data, objective_fun)

  % data
  has_constraint = data.has_constraint; 
  max_newton_iter = data.max_newton_iter;
  max_damp_iter = data.max_damp_iter;
  num_node = data.num_node;
  num_para = data.num_para; 
  gamma_d = data.gamma_d; 
  min_d = data.min_d; 
  if has_constraint
      u_old = data.u_old;
  else
      u_old = data.u_old(1:num_node);
  end
  tol = 1.0e2 * eps;
  fprintf('Function my_newton(): tol = %E ------ \n', tol); 

  objective = zeros(max_newton_iter+1, 1);
  norm_jacobian  = zeros(max_newton_iter+1, 1);
  norm_du = zeros(max_newton_iter+1, 1);
  norm_dd = zeros(max_newton_iter+1, 1);
  
  if has_constraint
        % data.Cv_bar and data.d0 are only used inside the my_newton function, 
        % but not outside. 
        data.Cv_bar = zeros(num_node, num_para);
        last_string = 'Diffusion Coefficients';
  else
        last_string = 'x';
  end
  % Newton iterations
  % Objective Function, Norm of Jacobian, Damping Parameter
  % Residual, Diffusion Coefficients
  fprintf('Obj Func \t Norm Jac \t Damping \t ');
  fprintf('Residual \t %s\n', last_string);
  %
  norm_du(1) = norm(u_old(1:num_node) - data.u2);
  norm_dd(1) = 0; 
  d_old = data.d0; 
  for i = 1: max_newton_iter
        if i == 1
            data.x = u_old; 
            data.d = d_old; 
            % The jacobian is calculated first to update the matrices. 
            data = objective_fun(data, 1); 
            jacobian_old = data.y;
            norm_jacobian(i) = norm(jacobian_old);
            data = objective_fun(data, 0);
            objective(i) = data.y;
        end
        
        % Calculate Hessian
        data = objective_fun(data, 2); 
        hessian_old = data.y;

        if norm_jacobian(i) < tol
%            u_new = u_old; 
%            d_new = d_old; 
           break;
        end

        switch data.cell_name
            case 'general_test'
                du = - hessian_old \ jacobian_old; 
            otherwise
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
            du = [data.delta_u; data.delta_v]; 
            delta_d = data.delta_d; 
        end
        %

        if ~max_damp_iter % no damping
            ss = 1; 
            u_new = u_old + ss * du;
            d_new = d_old + ss * delta_d; 
%             data.x = u_new;
%             data.d = d_new;
%             data = objective_fun(data, 1); 
%             jacobian_new = data.y;
        else
            % Damping
            alpha0 = norm(hessian_old * [du; delta_d] + jacobian_old) / norm_jacobian(i);
            delta = (1 - alpha0) / 2; % set delta: 0 < delta < 1 - alpha0
            % KK = data.KK; 
            KK = 0;
            for damp_iter = 1:max_damp_iter
                ss = 1.0/(1+KK * norm_jacobian(i));
                %%%%%%%%%%%%%%%%%%%%%
                u_new = u_old + ss * du;
                data.x = u_new;
                d_new = d_old + ss * delta_d; 
                data.d = d_new;
                data = objective_fun(data, 1);
                jacobian_new = data.y;
                %%%%%%%%%%%%%%%%%%%%%
                if (1 - norm(jacobian_new)/norm_jacobian(i))/ss < delta
                  if KK == 0
                    KK = 1;
                  else
                    KK = 10 * KK;
                  end
                else
                    data.KK = KK / 10; % Kathy: should this be KK, not data.KK? 
                    break;
                end
            end % for damp_iter = 1:max_damp_iter
            % fprintf('# of damping steps at step %d: %d \n', i, damp_iter);
        end % if ~max_damp_iter % no damping

        if has_constraint
            % constraint dd >= min_d
            % This step is necessary in some cases. 
            % This bound is only optionally triggered. Need to print something
            % here. 
            d_new = max(d_new, min_d);
%             uu = u_new(1:num_node);
%             dd = d_new; 
            print_para = min(data.num_para, 5); % print the first 5 parameters
        else
%             uu = u_new(1:num_node); 
            d_new = 0; 
        end
        % calculate objective function and Jacobian
        data.x = u_new;
        data.d = d_new; 
        data = objective_fun(data, 1); % Jacobian needs to calculate first to update matrices
        jacobian_new = data.y; 
        norm_jacobian(i+1) = norm(jacobian_new);
        data = objective_fun(data, 0);
        objective(i+1) = data.y;
        
        norm_du(i+1) = norm(u_new(1:num_node) - data.u2);
        norm_dd(i+1) = norm((d_new - d_old)./data.d0);

      % From left to right: Objective Function, Norm of Jacobian, Damping Parameter
      % Residual, Diffusion Coefficients

        fprintf('%8.3e \t %10.3e \t %5.2e \t %5.2e \t ', ...
            objective(i), norm_jacobian(i), ss, norm_du(i));
        if has_constraint % print diff_coef
            fprintf('%s \n', num2str(d_old(1:print_para)'));
            % convergence test
            tmp = (d_new - d_old)./d_old; 
            criteria = (norm_jacobian(i+1) < sqrt(tol) && max(abs(tmp)) < 0.01);
            clear tmp; 
        else
            fprintf('%s \n', num2str(u_old'));
            criteria = (norm_jacobian(i+1) < tol); 
        end

        if criteria % True: (norm_jacobian(i+1) < 1e-06 && max(abs(tmp)) < 0.01)
            fprintf('%8.3e \t %10.3e \t %5.2e \t %5.2e \t ', ...
                objective(i+1), norm_jacobian(i+1), 0, norm_du(i+1));
            if has_constraint % print diff_coef
                fprintf('%s \n', num2str(d_old(1:print_para)'));
            else
                fprintf('%s \n', num2str(u_new'));
            end
            break; 
        else
            clear u_old d_old jacobian_old;
            u_old = u_new;
            d_old = d_new; 
            jacobian_old = jacobian_new; 
        end        
  end % for i = 1 : max_newton_iter

  outer_iter = data.outer_iter;
  num_newton_iter = i+1;
  if has_constraint
    data.u_new = u_new;
    data.d_old = d_old; 
    data.d_new = d_new; 
  else
      vv = zeros(num_node, 1);
      data.u_new = [u_new; vv];
      clear vv; 
  end
  data.norm_du{outer_iter} = norm_du(1:num_newton_iter);
  data.norm_dd{outer_iter} = norm_dd(1:num_newton_iter)*sqrt(num_node*gamma_d);

  fprintf('Number of newton steps: %d \n', i);
  fprintf('######################## \n\n');
  
return

