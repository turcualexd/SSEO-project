clear
close all
clc

string_short = 848;
cell_series_short = 13;
parallel_short = 64;

string_medium = 369;
cell_series_medium = 14;
parallel_medium = 40;

string_long = 114;
cell_series_long = 22;

string_total = string_long + string_medium + string_short; % --- valore coerente ---

cell_weight = 84; % mg/cm2
cells_total = 18698; % --- juno press kit ---- 
cells_total_calc = string_short * cell_series_short + string_medium * cell_series_medium + string_long * cell_series_long; 
cells_weight_tot = 26.6 * 1e-4 * cells_total * cell_weight * 1e-3;
% --- valore calcolato coerente --- 

% ---- Vmp at 5.44 AU ---- 
Vmp_short_544 = 36; 
Vmp_medium_544 = 39;
Vmp_long_544 = 62;

Vmp_cell_short_544 = Vmp_short_544/cell_series_short; % 2.7692
Vmp_cell_medium_544 = Vmp_medium_544/cell_series_medium; % 2.7857
Vmp_cell_long_544 = Vmp_long_544/cell_series_long; % 2.8182

% ---- Vmp at 1 AU ---- 
% Vmp_short_1 = 36; 
% Vmp_medium_1 = 39;
% Vmp_long_1 = 62;
% 
% Vmp_cell_short_1 = Vmp_short/cell_series_short; % 2.7692
% Vmp_cell_medium_1 = Vmp_medium/cell_series_medium; % 2.7857
% Vmp_cell_long_1 = Vmp_long/cell_series_long; % 2.8182

total_area = 8 * 2.023 * 2.7 + 3 * 2.023 * 2.374;

area_cell = total_area/cells_total * 1e4; % 31.0753 cm2 valore coerente con le celle trovate






