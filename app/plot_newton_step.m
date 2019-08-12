function plot_newton_step(y, label_str)
my_figure('font_size', 24, 'line_width', 3);  
num_ns = length(y); % number of Newton steps
x = (0:num_ns-1)';
plot(x, y, '-o', 'LineWidth', 2); 
% strValue = strtrim(cellstr(num2str([x y'], '(%d, %4.3e)')));
% text(x,y, strValue, 'VerticalAlignment', 'bottom', 'FontSize', 12);
temp = axis;
axis([-1 num_ns temp(3) temp(4)]);
xticks(-1:1:num_ns);
xlabel('Newton Step'); ylabel(label_str); 
return