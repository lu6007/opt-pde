function data = hessian_sub_mat(data, u, v)

    cell_name = data.cell_name;
    num_para = data.num_para;
    d0 = data.d0;
    A_sub = data.A_sub;
    M_sub = data.M_sub;
    gamma = data.gamma_d; 

  % different model problem should be dealt with differently
  coef_1 = {'layered_diffusion','spot_diffusion', 'mem17', 'test_cont','WTH3K9', ...
      'layered_diffusion_general'};
  % coef_3 = {'tensor_diffusion','tensor_cross','tensor_cross_2'};

  num_row = size(A_sub{1},1);
  Cv = zeros(num_row, num_para);
  Cu = zeros(num_row, num_para);
  DG = zeros(num_para, 1);
  if sum(ismember(coef_1, cell_name)) > 0

    % format linear system

    for i = 1 : num_para
      Cv(:, i) = A_sub{i} * u;
      Cu(:, i) = A_sub{i} * v;
      if abs(d0(i)) > eps
        DG(i)    = sum(M_sub{i}(:))/(d0(i)^2); 
      else
        DG(i)    = sum(M_sub{i}(:)); 
      end
    end
    G = gamma * diag(DG);

  else
    switch cell_name
      case 'tensor_diffusion'

        % formulate the linear system
        Cu(:, 1) = A_sub{1} * v;
        Cv(:, 1) = A_sub{1} * u;
        Cu(:,2) = (A_sub{2}+A_sub{4})*v;
        Cv(:,2) = (A_sub{2}+A_sub{4})*u; 
        Cu(:,3) = (A_sub{3}+A_sub{5})*v;
        Cv(:,3) = (A_sub{3}+A_sub{5})*u;

        DG(1) = gamma*sum(M_sub{1}(:))/(d0(1)^2);
        DG(2) = gamma*0.5*(sum(M_sub{2}(:))+sum(M_sub{3}(:)))/(d0(2)^2);
        DG(3) = gamma*0.5*(sum(M_sub{2}(:))+sum(M_sub{3}(:)))/(d0(3)^2); 

        G = diag(DG);

      case 'tensor_cross_2'

        % formulate the linear system
        Cu(:, 1) = A_sub{1} * v;
        Cv(:, 1) = A_sub{1} * u;
        Cu(:,2) = (A_sub{2} +A_sub{4}+A_sub{6})*v;
        Cv(:,2) = (A_sub{2}+A_sub{4}+A_sub{6})*u; 
        Cu(:,3) = (A_sub{3}+A_sub{5})*v;
        Cv(:,3) = (A_sub{3}+A_sub{5})*u;

        DG(1) = gamma*sum(M_sub{1}(:))/(d0(1)^2);
        DG(2) = gamma*0.5*(sum(M_sub{2}(:))+sum(M_sub{3}(:))+sum(M_sub{4}(:)))/(d0(2)^2);
        DG(3) = gamma*0.5*(sum(M_sub{2}(:))+sum(M_sub{3}(:)))/(d0(3)^2);

        G = diag(DG);
    end % switch cell_name
  end % if sum(ismember(coef_1, cell_name)) > 0
  
  data.Cu = Cu;
  data.Cv = Cv;
  data.G = G; 
  
  clear Cu Cv G;
  
  return
