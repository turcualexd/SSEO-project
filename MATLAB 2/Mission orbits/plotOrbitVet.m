function plotOrbitVet(rvet, n, color, tipo)
if nargin == 4
    plot3(rvet(1, :), rvet(2, :), rvet(3, :), tipo, 'LineWidth', n, 'Color', color);
elseif nargin == 3
    plot3(rvet(1, :), rvet(2, :), rvet(3, :), 'LineWidth', n, 'Color', color);
elseif nargin == 2
    plot3(rvet(1, :), rvet(2, :), rvet(3, :), 'LineWidth', n);
else
    plot3(rvet(1, :), rvet(2, :), rvet(3, :));
end
