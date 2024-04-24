clear
clc
close all

juno = load("juno_table_jupiter.txt");


for i = 1 : size(R_juno_jupiter, 2)
    B(i) = 2*4.3e-4 * 71372^3 / norm(R_juno_jupiter(:,i))^3;
end

Mm = B * 20;

R_juno_sun = loadData("horizons_results_sun.txt");


for i = 1 : size(R_juno_sun, 2)
    R_js =  norm(R_juno_sun(:,i));
end

n = 2 * 2 * pi / 60; % velocità attorno a giove