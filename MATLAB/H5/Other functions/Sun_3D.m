function Sun_3D(n)

%---------------------------Plot Sun in 3D plot----------------------------
%
% This function plots a sphere with Sun texture in a pre-existing 3D space
% at origin.
%
%--------------------------------------------------------------------------
%
% INPUT:
% n                 [1]: scale factor for Sun radius (default 1)
%
%--------------------------------------------------------------------------
%
% AUTHOR: Alex Cristian Turcu
%
%--------------------------------------------------------------------------

%% Default Input

if nargin < 1
    n = 1;
end
Rt = astroConstants(3) * n / astroConstants(2);

image = 'Sun.jpg';

%% Figure

% Create the figure
hold on;
grid on;

% Set the axes scale equal
axis equal;

% Put the axes labels
xlabel('X [AU]', 'FontSize', 15);
ylabel('Y [AU]', 'FontSize', 15);
zlabel('Z [AU]', 'FontSize', 15);

% Set initial view
view(120,30);

% Define the number of panels to be used to model the sphere 
npanels = 180;  

% Create a 3D meshgrid of the sphere points using the ellipsoid function
[x, y, z] = ellipsoid(0, 0, 0, Rt, Rt, Rt, npanels);

% Create the globe with the surf function
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');

%% Texturemap the globe

cdata = imread(image);

% Set the transparency of the globe: 1 = opaque, 0 = invisible
alpha = 1; 

% Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
% specify the image data using the 'CData' property with the data loaded 
% from the image. Finally, set the transparency and remove the edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

end