function drawPlanet(planet, position, scale, f)

if nargin == 3
    f = 0;
elseif nargin == 2
    scale = 1;
    f = 0;
end
if planet > 11 && planet < 100
    planet = -1;
end

% Directory of textures folder
img = 'Planet textures';

% Select planet
switch planet
    case 10
        R = astroConstants(3);
        img = strcat(img, '\Sun.jpg');
    case 1
        R = astroConstants(21);
        img = strcat(img, '\Mercury.jpg');
    case 2
        R = astroConstants(22);
        img = strcat(img, '\Venus.jpg');
    case 3
        R = astroConstants(23);
        img = strcat(img, '\Earth.jpg');
    case 4
        R = astroConstants(24);
        img = strcat(img, '\Mars.jpg');
    case 5
        R = astroConstants(25);
        img = strcat(img, '\Jupiter.jpg');
    case 6
        R = astroConstants(26);
        img = strcat(img, '\Saturn.jpg');
    case 7
        R = astroConstants(27);
        img = strcat(img, '\Uranus.jpg');
    case 8
        R = astroConstants(28);
        img = strcat(img, '\Neptune.jpg');
    case 11 
        R = astroConstants(30);
        img = strcat(img, '\Moon.jpg');
    case -1
        R = astroConstants(23);
        img = strcat(img, '\Asteroid.jpg');
    otherwise
        error("Wrong planet ID")
end
R = scale*R;

% Set figure
if f > 0
    background_plot = 'w';
    figure('Color', background_plot);
    hold on;
    grid minor;
    axis equal;
    xlabel('X [km]');
    ylabel('Y [km]');
    zlabel('Z [km]');
    view(120,30);
end

% Define wireframe
npanels = 180;  
[x, y, z] = ellipsoid(position(1), position(2), -position(3), R, R, R, npanels);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');

% Texture the planet
cdata = imread(img);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');

end