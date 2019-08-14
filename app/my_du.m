function data = my_du(data, jacobian_old)

  num_node = data.num_node; 

  % assemble A and M
  A = data.A; M = data.M; 
  Cu = data.Cu; Cv = data.Cv; G = data.G;
  Cv_bar = data.Cv_bar; 
  Ju = jacobian_old(1:num_node); 
  Jv = jacobian_old(num_node+1: 2*num_node);
  Jd = jacobian_old(2*num_node+1: end);

  Cv_hat = Cv_bar; clear Cv_bar;
  %
  Jv_bar = A \ Jv;
  W = A \ (Cv - A * Cv_hat);
  Cv_bar = Cv_hat + W; 
  Ju_bar = Ju - M * Jv_bar;
  Cu_bar = Cu - M * Cv_bar;
  G_bar = G - Cu' * Cv_bar - Cv_bar' * Cu_bar;
  % delta_d = - pinv(G_bar) * (Jd - Cu' * Jv_bar - Cv_bar' * Ju_bar);
  
  %
  data.Cv_bar = Cv_bar;
  data.delta_d = - G_bar \ (Jd - Cu' * Jv_bar - Cv_bar' * Ju_bar);
  data.delta_u = - (Jv_bar + Cv_bar * data.delta_d);
  data.delta_v = - A \ (Ju_bar + Cu_bar * data.delta_d);
end
