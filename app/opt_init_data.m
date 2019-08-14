function data = opt_init_data(name)
switch name
    case 'spot'
        data.max_newton_iter = 10;
        data.init_u_tag = 1; % u0 = u1 or u2
        data.init_d = 10;
        data.gamma = 1.0e-5;
        data.enable_normalize = 0;
        data.enable_damp_newton = 0;
        data.max_damp_step = 10;
        data.max_outer_iter = 10;
        data.fun_idx = 7;
    case 'layer'
        data.max_newton_iter = 20;
        data.init_u_tag = 1; % u0 = u1 or u2
        data.init_d = 10;
        data.gamma = 1.0e-5;
        data.enable_normalize = 1;
        data.enable_damp_newton = 0;
        data.max_damp_step = 10;
        data.max_outer_iter = 10;
        data.fun_idx = 7;
    case 'tensor'
        data.max_newton_iter = 10;
        data.init_u_tag = 1; % u0 = u1 or u2
        data.init_d = 30; % 30 not convergent with or without damping
        data.gamma = 1.0e-5; % 1.0, 0.1
        data.enable_normalize = 1;
        data.enable_damp_newton = 1;
        data.max_damp_step = 10;
        data.max_outer_iter = 1;
        data.fun_idx = 7;
    case 'tensor_cross_2'
        data.max_newton_iter = 20;
        data.init_u_tag = 1; % u0 = u1 or u2
        data.init_d = 10; % 30 not convergent with or without damping
        data.gamma = 1e-5; % 1.0, 0.1
        data.enable_normalize = 0; % for 4, normalizaton and damping not work
        data.enable_damp_newton = 0;
        data.max_damp_step = 10;
        data.max_outer_iter = 1;
        data.fun_idx = 7;
    case 'mem17'
        data.max_newton_iter = 10;
        data.init_u_tag = 2; % u0 = u1 or u2
        data.init_d = 60; 
        data.gamma = 1e-5; % 1.0, 0.1
        data.enable_normalize = 0;
        data.enable_damp_newton = 0;
        data.max_damp_step = 10;
        data.max_outer_iter = 5;
        data.fun_idx = 7;
    case {1, 2, 3, 4, 5}
        data.max_newton_iter = 20;
        data.init_u_tag = 2; % u0 = u1 or u2
        data.init_d = 10; % 30 not convergent with or without damping
        data.gamma = 1e-5; % 1.0, 0.1
        data.enable_normalize = 0;
        data.enable_damp_newton = 0;
        data.max_damp_step = 10;
        data.max_outer_iter = 1;
        data.fun_idx = name;
    case 6
        data.max_newton_iter = 20;
        data.init_u_tag = 2; % u0 = u1 or u2
        data.init_d = 10; 
        data.gamma = 1e-5; 
        data.enable_normalize = 0;
        data.enable_damp_newton = 1;
        data.max_damp_step = 10;
        data.max_outer_iter = 1;
        data.fun_idx = name;        
    case 'general'
        data.max_newton_iter = 5;
        data.init_u_tag = 2;
        data.init_d = 10;
        data.gamma = 1e-5;
        data.enable_normalize = 0;
        data.enable_damp_newton = 0;
        data.max_damp_step = 10;
        data.max_outer_iter = 1;
        data.fun_idx = 7;
%     case {'layered_diffusion_general', 9}
%         data.max_newton_iter = 10;
%         data.init_u_tag = 2;
%         data.init_d = 10;
%         data.gamma = 1e-5;
%         data.enable_normalize = 1;
%         data.enable_damp_newton = 1;
%         data.max_damp_step = 10;
%         data.max_outer_iter = 10;
%         data.fun_idx = 8;
end
data.load_data_function = @load_data;
data.normalize_data_function = @normalize_data;
return;

function data = load_data(name)
  pa = '../data/';
  % pa = '~/Desktop/opt_pde_0716/data/';
  switch name
    case 'spot'
      data_file = strcat(pa, 'spot_diffusion.mat');
      data = load(data_file);
      data.tri(4, data.tri(4, :) == 3) = 2;
      data.tri(4, data.tri(4, :) == 4) = 2;
      data.tri(4, data.tri(4, :) == 5) = 2;
      data.cell_name = 'spot_diffusion';
      data.num_para = size(unique(data.tri(4, :)), 2);
      data.coef_tag = 1;
      data.fem = 1;
      data.default_solver = 0;
    case 'layer'
      data_file = strcat(pa, 'layered_diffusion.mat');
      data = load(data_file);
      data.cell_name = 'layered_diffusion';
      data.num_para = size(unique(data.tri(4, :)), 2);
      data.coef_tag = 1;
      data.fem = 1;
      data.default_solver = 0;
    case 'tensor'
      data_file = strcat(pa, 'tensor_diffusion.mat');
      data = load(data_file);
      data.cell_name = 'tensor_diffusion';
      data.num_para = size(unique(data.tri(4, :)), 2);
      data.coef_tag = 2;
      data.fem = 1;
      data.default_solver = 0;
      % data.KK = 0;
    case 'tensor_cross_2'
      data_file = strcat(pa, 'tensor_cross_2.mat');
      data = load(data_file);
      data.tri(4, data.tri(4, :) == 3) = 2;
      data.tri(4, data.tri(4, :) == 4) = 3;
      data.tri(4, data.tri(4, :) == 5) = 3;
      data.tri(4, data.tri(4, :) == 6) = 4;
      data.cell_name = 'tensor_cross_2';
      data.num_para = 3;
      data.coef_tag = 2;
      data.fem = 1;
      data.default_solver = 0;
    case 'mem17'
      data_file = strcat(pa, 'mem17.mat');
      data = load(data_file);
      data.cell_name = 'mem17';
      data.u1 = data.u(:, 1);
      data.u2 = data.u(:, 2);
%       data.p = data.p_image;
%      data.p_image = data.p;
      data.num_para = 1;
      data.coef_tag = 1;
      data.fem = 1;
      data.default_solver = 0;
      data.plot_surf = 1;
    case {1, 2, 3, 4, 5, 6}
      data.cell_name = 'general_test';
      data.num_node = 1;
      data.num_para = 0;
      % data.fun_idx = name;
      data.u0 = 1;
      data.u2 = 0;
      data.coef_tag = 1;
      data.fem = 0;
      data.default_solver = 1;
    case 'general'
      data_file = strcat(pa, 'layered_diffusion_general_5_refined 2.mat');
      data = load(data_file);
      data.cell_name = 'layered_diffusion_general';
      for i = 1 : size(data.tri, 2)
        data.tri(4, i) = i;
      end
      data.num_para = size(unique(data.tri(4, :)), 2);
      data.coef_tag = 1;
      data.p = data.p_image;
      data.fem = 1;
      data.default_solver = 0;
      data.plot_surf = 1;
%     case {9, 'layered_diffusion_general'}
%       data_file = strcat(pa, 'layered_diffusion_general_5_refined.mat');
%       data = load(data_file);
%       data.cell_name = 'layered_diffusion_general_h1';
%       for i = 1 : size(data.tri, 2)
%         data.tri(4, i) = i;
%       end
%       % data.num_para = size(unique(data.tri(4, :)), 2);
%       data.num_para = size(data.p_image, 2);
%       data.coef_tag = 1;
%       data.p = data.p_image;
%       data.fem = 1;
%       data.default_solver = 0;
%       data.plot_surf = 1;
  end

  if data.fem == 1
    % change tri, edge, p, and p_image to column vectors
    tri = (data.tri)';
    edge = (data.edge)';
    p = (data.p)';
%    p_image = (data.p_image)';
    data.num_node = size(data.p, 2);
    data = rmfield(data, {'tri', 'edge','p', 'p_image'});
    data.tri = tri;
    data.edge = edge;
    data.p = p;
%     data.p_image = p_image;
%    data.p_image = p; %%% Kathy 08/01/2018 This is weird
    clear tri edge p p_image;

    data.path = pa;
%     data.num_para = size(unique(data.tri(:,4)), 1);
    if strcmp(name, 'tensor_cross_2')
        data.num_para = 3;
    end

    % Remove some unnecessary fields
    data = rmfield(data, 'num_nodes');
%     if isfield(data, 'coef_tag') % needed for name == 5
%         data = rmfield(data, 'coef_tag');
%     end
    data = rmfield(data, 'u');
  end

return;

function data = normalize_data(data, enable_normalize)
  if enable_normalize
      % --------------- normalization ---------------
      % zmax = max(max(data.u1),max(data.u2));
      % zmin = min(min(data.u1),min(data.u2));
      abs_u_diff = abs(data.u2-data.u1);
      umax = max(abs_u_diff);
      umin = min(abs_u_diff);
      xmax = max(data.p(:,1));
      xmin = min(data.p(:,1));
      ymax = max(data.p(:,2));
      ymin = min(data.p(:,2));
      xmid  = (xmax+xmin)/2.0;
      ymid  = (ymax+ymin)/2.0;
      xy_scale = 1.0 / max(xmax-xmin,ymax-ymin);
      u_scale = 1.0 / (umax-umin);

      p1 = xy_scale * (data.p(:, 1) - xmid);
      p2 = xy_scale * (data.p(:, 2) - ymid);
  %     data.u1 = (data.u1 - umin)*u_scale;
  %     data.u2 = (data.u2 - umin)*u_scale;
      data.u1 = data.u1 * u_scale;
      data.u2 = data.u2 * u_scale;
      data.p = [p1, p2];
      data.dt = data.dt * (xy_scale^2);
      data.gamma_d = data.gamma_d * xy_scale;
      % init_d = init_d * xy_scale;

      % --------------- end ---------------
  else
      xy_scale = 1;
  end
  data.xy_scale = xy_scale;

return;
