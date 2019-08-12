obj_tmp = cell2mat(data.objective');
s_tmp = cell2mat(data.s_hist');
d_tmp = cell2mat(data.d_hist) / xy_scale;
% d_tmp = d_tmp' ./ [29, 5];
% d_tmp = d_tmp' ./ [5, 20, 30, 1];
d_tmp = d_tmp' ./ [5, 100, 5];
xtick_step = 1;


figure;
plot(obj_tmp, '.-', 'MarkerSize', 14); 
% title('objective, normalized', 'FontSize', 16);
title('objective, normalized, damped', 'FontSize', 16);
% title('objective, normalized', 'FontSize', 16);
xticks(1 : xtick_step : size(obj_tmp, 1))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',20)
grid on

figure;
plot(s_tmp, '.-', 'MarkerSize', 14); 
% title('objective, normalized', 'FontSize', 16);
title('damping parameter, normalized, damped', 'FontSize', 16);
% title('objective, normalized', 'FontSize', 16);
xticks(1 : xtick_step : size(obj_tmp, 1))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',20)
grid on

figure; 
plot(d_tmp, '.-', 'MarkerSize', 14); 
% title('diffusion coefficient, normalized', 'FontSize', 16);
title('diffusion coefficient, normalized, damped', 'FontSize', 16);
% title('diffusion coefficient, normalized', 'FontSize', 16);
xticks(1 : xtick_step :size(obj_tmp, 1))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',20)
% legend('29', '5')
% legend('5', '20', '30', '1')
% legend('5', '100', '5')
grid on
