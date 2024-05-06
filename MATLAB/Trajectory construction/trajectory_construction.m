clear; clc; close all;

%% 
% converte linee multiple di file .txt
F = fopen("EARTH_horizon_long.txt", "r");
G = fopen("EARTH_long.txt", "w");
if F ~= -1 && G ~= -1
    while ~feof(F)
        line = fgetl(F);
        year = str2double(line(26:29));
        switch line(31:33)
            case "Jan"
                month = 1;
            case "Feb"
                month = 2;
            case "Mar"
                month = 3;
            case "Apr"
                month = 4;
            case "May"
                month = 5;
            case "Jun"
                month = 6;
            case "Jul"
                month = 7;
            case "Aug"
                month = 8;
            case "Sep"
                month = 9;
            case "Oct"
                month = 10;
            case "Nov"
                month = 11;
            case "Dec"
                month = 12;
        end
        day = str2double(line(35:36));
        hrs = str2double(line(38:39));
        min = str2double(line(41:42));
        sec = str2double(line(44:50));
        date = date2mjd2000([year month day hrs min sec]);
        line = fgetl(F);
        e = str2double(line(6:26));
        i = str2double(line(58:78));
        line = fgetl(F);
        OM = str2double(line(6:26));
        om = str2double(line(32:52));
        line = fgetl(F);
        theta = str2double(line(58:78));
        line = fgetl(F);
        a = str2double(line(6:26));
        fprintf(G, "%f %f %f %f %f %f %f\n", date, a, e, i, OM, om, theta);
    end
end

fclose(G);
fclose(F);

%%

ch = 16;
lw = 2;
sun = 30;
earth = 600;
assi = 16;

el = load("EARTH_long.txt");
el = [el(:, 2:3) deg2rad(el(:, 4:7))];
mu = astroConstants(4);
rvet = zeros(3, size(el, 1));
for i = 1 : length(rvet)
    rvet(:, i) = kep2car(el(i, :), mu);
end
% plot for 30days
rvet = rvet(:, 1:60);
plotOrbitVet(rvet, lw, "#41bf13")
hold on
grid minor
axis equal
drawPlanet(10, [0 0 0], sun, 0)
drawPlanet(3, rvet(:, 1), earth, 0);
drawPlanet(3, rvet(:, end), earth, 0);


%legend("Juno trajectory", "", "", "", "", "Earth orbit", "Mars orbit", "Jupiter orbit", 'FontSize', ch)
xlabel("[km]", 'Interpreter', 'latex', 'FontSize', assi)
ylabel("[km]", 'Interpreter', 'latex', 'FontSize', assi)