clear
clc
close all

R = loadData("horizons_results_jupiter.txt");

for i = 1 : size(R, 2)
    B(i) = 2*4.3e-4 * 71372^3 / norm(R(:,i))^3;
end

Mm = B * 20;