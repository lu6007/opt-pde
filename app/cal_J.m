function [Ju, Jv, Jd] = cal_J(data, d, u, v)
    cell_name = data.cell_name;
    num_para = data.num_para;
    d0 = data.d_old;
    u2 = data.u2;
    u1 = data.u1;
    gamma = data.gamma_d;
    dt = data.dt;
    A = data.A;
    M = data.M;
    A_sub = data.A_sub;
    M_sub = data.M_sub;

  % different model problem should be dealt with differently
  coef_1 = {'layered_diffusion','spot_diffusion', 'mem17', 'test_cont','WTH3K9', ...
      'layered_diffusion_general'};
  % coef_3 = {'tensor_diffusion','tensor_cross','tensor_cross_2'};

  Jd = zeros(num_para, 1);
  if sum(ismember(coef_1, cell_name)) > 0

    % evaluate J for current [u; v; d]
    Ju = M * (u - u2) + A * v;
    Jv = A * u - M * u1 / dt;
    % Jv = A * u - M * (u - u1)/dt; 

    dd0 = d - d0;
    for i = 1 : num_para      
      if abs(d0(i)) > eps
        Jd(i, 1) = gamma * sum(M_sub{i}(:)) * dd0(i) / (d0(i)^2) + ...
          u' * A_sub{i} * v;
      else % abs(d0(i)) <= eps
        Jd(i, 1) = gamma * sum(M_sub{i}(:)) * dd0(i) + u' * A_sub{i} * v;
      end
    end

  else

    Ju = M * (u - u2) + A * v;
    Jv = A * u - M * u1 / dt;
    % Jv = A*u + M * (u-u1) / dt;
    dd0 = d - d0;
    switch cell_name
      case 'tensor_diffusion'

        Jd(1, 1) = gamma * sum(M_sub{1}(:)) * dd0(1) / (d0(1)^2) + u' * A_sub{1} * v;
        Jd(2, 1) = gamma * sum(M_sub{2}(:)) * dd0(2) / (d0(2)^2) + u' * (A_sub{2}+A_sub{4}) * v;
        Jd(3, 1) = gamma * sum(M_sub{3}(:)) * dd0(3) / (d0(3)^2) + u' * (A_sub{3}+A_sub{5}) * v;

      case 'tensor_cross_2'

        Jd(1, 1) = gamma * sum(M_sub{1}(:)) * dd0(1) / (d0(1)^2) + u' * A_sub{1} * v;
        sumM2 = sum(M_sub{2}(:));
        sumM3 = sum(M_sub{3}(:));
        sumM23_mean = 0.5 *(sumM2+sumM3);
        Jd(2, 1) = gamma * sumM23_mean * dd0(2) / (d0(2)^2) + u' * (A_sub{2}+A_sub{4}+A_sub{6}) * v;
        Jd(3, 1) = gamma * sumM23_mean * dd0(3) / (d0(3)^2) + u' * (A_sub{3}+A_sub{5}) * v;

    end
  end
