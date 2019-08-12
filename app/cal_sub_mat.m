% Calculate the sumatrices and update diffusion coefficients
% with different cases
function data = cal_sub_mat(data, d)

    coef_1 = {'layered_diffusion','spot_diffusion', 'mem17', 'test_cont', 'WTH3K9', ...
        'layered_diffusion_general'};
    coef_3 = {'tensor_diffusion','tensor_cross','tensor_cross_2'};
    cell_name= data.cell_name;
    num_para = data.num_para; 
    tri = data.tri;
    p = data.p;
    
    % Stiffness and mass matrix for different subregion
    A_sub = cell(num_para, 1);
    M_sub = cell(num_para, 1);
    
  if sum(ismember(coef_1, cell_name)) > 0
    diff_coef = ones(1, size(tri, 2));
    for i = 1 : num_para
      idx_tmp = (tri(4, :) == i);
      tri_tmp = tri(:, idx_tmp);
      
      % [A_sub{i}, M_sub{i}, ~] = assema(p, tri_tmp, diff_coef(idx_tmp), 1, 0);
      [A_sub{i}, M_sub{i}] = assemble_matrix(p, tri_tmp, 'diff_coef', 1);
      diff_coef(idx_tmp) = d(i); 
      clear indx_tmp tri_tmp; 
    end

  elseif sum(ismember(coef_3, cell_name)) > 0

        % load orientation data
        orientation_file = strcat(data.path, cell_name, '/', 'orientation.mat');
        ori = load(orientation_file);
        x   = ori.x;
        y   = ori.y;
        clear ori;

        num_line = numel(x)/2;
        TT = cell(num_line, 1); 
        for i = 1 : num_line
              a = x(2*i) - x(2*i-1);
              b = y(2*i) - y(2*i-1);
              % Normalizing
              temp = sqrt(a^2 + b^2);
              a = a / temp;
              b = b / temp;

              % Calculate the diffusion tensor with the transforming matrix T
              TT{i} = [ a b; ...
                       -b a]';
              clear temp a b;
        end

    switch cell_name
      case 'tensor_diffusion'
        diff_coef = ones(3, size(tri, 2));
        diff_const = d;

        num_D = num_para+ num_line;
        D_tmp = cell(num_D,1);
        D_tmp{1} = [1, 0; 0, 1];
        for i = 1:num_line
          T = TT{i};
          D_tmp{2+2*(i-1)} = T'*[1 0; 0 0]*T;
          D_tmp{3+2*(i-1)} = T'*[0 0; 0 1]*T;
        end

        D = cell(num_D,1);
        for i = 1:num_D
            D{i} = [D_tmp{i}(1,1); D_tmp{i}(1,2); D_tmp{i}(2,2)];
        end
        clear D_tmp;

		% mass and stiffness matrix for rest of the cell region
        idx_tmp = (tri(4, :) == 1);
        tri_tmp = tri(:, idx_tmp);
        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{1});
%        [A_sub{1}, M_sub{1}, ~] = assema(p, tri_tmp, 1, 1, 0);
        [A_sub{1}, M_sub{1}] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:,idx_tmp));

        % mass and stiffness matrix for tensor 1 region, should have 2 sets of matrices
        idx_tmp = (tri(4, :) == 2);
        tri_tmp = tri(:, idx_tmp);
        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{2});
        [A_sub{2}, M_sub{2}] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));
        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{3});
        [A_sub{3}, ~] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));

        % ------- tensor 2 region

        idx_tmp = (tri(4, :) == 3);
        tri_tmp = tri(:, idx_tmp);
        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{4});
        [A_sub{4}, M_sub{3}] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));
        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{5});
        [A_sub{5}, ~] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));

		% ------- Update diffusion coefficients

        cartesian_diff_const = { ...
            diff_const(1)*D{1};...
            diff_const(2)*D{2}+diff_const(3)*D{3};...
            diff_const(2)*D{4}+diff_const(3)*D{5}};
        for i = 1:num_para
            idx_tmp = (tri(4,:)==i);
            diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, cartesian_diff_const{i});
        end

      case 'tensor_cross_2'
        diff_coef = ones(3, size(tri, 2));
        diff_const = d;

        num_D = num_para+ num_line;
        D_tmp = cell(num_D,1);
        D_tmp{1} = [1, 0; 0, 1];
        D_tmp{6} = D_tmp{1};
        for i = 1:num_line
          T = TT{i};
          D_tmp{2+2*(i-1)} = T'*[1 0; 0 0]*T;
          D_tmp{3+2*(i-1)} = T'*[0 0; 0 1]*T;
        end

        D = cell(num_D,1);
        for i = 1:num_D
            D{i} = [D_tmp{i}(1,1); D_tmp{i}(1,2); D_tmp{i}(2,2)];
        end
        clear D_tmp;

		% mass and stiffness matrix for rest of the cell region
        idx_tmp = (tri(4, :) == 1);
        tri_tmp = tri(:, idx_tmp);

        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{1});
%        [A_sub{1}, M_sub{1}, ~] = assema(p, tri_tmp, 1, 1, 0);
        [A_sub{1}, M_sub{1}] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:,idx_tmp));
        
        % mass and stiffness matrix for tensor 1 region, should have 3 sets of matrices
        idx_tmp = (tri(4, :) == 2);
        tri_tmp = tri(:, idx_tmp);
        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{2});
        [A_sub{2}, M_sub{2}] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));
        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{3});
        [A_sub{3}, ~] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));

        % ------- tensor 2 region

        idx_tmp = (tri(4, :) == 3);
        tri_tmp = tri(:, idx_tmp);
        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{4});
        [A_sub{4}, M_sub{3}] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));
        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{5});
        [A_sub{5}, ~] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:, idx_tmp));
        
        % --- tensor cross region: 4 
        idx_tmp = (tri(4, :) == 4);
        tri_tmp = tri(:, idx_tmp);
        diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, D{1});
        [A_sub{6}, M_sub{4}] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(:,idx_tmp));

		% ------- Update diffusion coefficients

        cartesian_diff_const = { ...
            diff_const(1)*D{1};...
            diff_const(2)*D{2}+diff_const(3)*D{3};...
            diff_const(2)*D{4}+diff_const(3)*D{5}; ...
            diff_const(2)*D{1}};
        for i = 1:4 % num_region
            idx_tmp = (tri(4,:)==i);
            diff_coef(:,idx_tmp) = expand_diff_coef(idx_tmp, cartesian_diff_const{i});
        end
    end % swtich cell_name

  else
    diff_coef = ones(1, size(tri, 2));
    for i = 1 : num_para
      idx_tmp = (tri(4, :) == i);
      tri_tmp = tri(:, idx_tmp);
      [A_tmp, M_tmp] = assemble_matrix(p, tri_tmp, 'diff_coef', diff_coef(idx_tmp));
      diff_coef(idx_tmp) = d(i);
      A_sub{i} = A_tmp;
      M_sub{i} = M_tmp;
    end
  end
  
  data.A_sub = A_sub;
  data.M_sub = M_sub;
  data.diff_coef = diff_coef;
  
  clear A_sub M_sub diff_coef; 
  return
  
  function diff_coef = expand_diff_coef(index, tensor)
      % Copy a 3x1 tensor vector to a matrix of diffusion cefficients
      num_col = sum(index);
      diff_coef = zeros(3, num_col);
      for i = 1:3
          diff_coef(i,:) = tensor(i);
      end
  return
        
