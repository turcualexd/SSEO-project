function [t, r, v] = plotOrbit(r0, v0, varargin)

% plotOrbit - Plot an orbit in a 3D space.
%
% PROTOTYPE:
%   [t, r, v] = plotOrbit(r0, v0, deltaT, n, options)
%
% DESCRIPTION:
%   Function to propagate an orbit through the integration of the equation
%   of motion, given as input the initial radius and velocity.
%
% INPUT:
%   r0              [3x1]       Initial position vector             [km]
%   v0              [3x1]       Initial velocity vector             [km/s]
%   deltaT (opt.)   [1x1]       Time interval (negative values for
%                               the number of full revolutions)
%                               (Default: -1)                       [s]
%   n (opt.)        [1x1]       Size of the problem (vectors)
%                               (Default: 10000)                    [-]
%
% NAME-VALUE ARGUMENTS:
%   'mu'            [1x1]       Standard gravitational parameter
%                               (Default: Earth -> 3.98600433e+5)   [km^3/s^2]
%   'color'         [char]      Color of the orbit plot
%                               (Default: 'b')
%
% OUTPUT:
%   (plot3)         [3D line]   Plots the orbit on an existing figure
%   t               [nx1]       Times vector                        [s]
%   r               [nx3]       Position vectors (for each time)    [km]
%   v               [nx3]       Velocity vectors (for each time)    [km/s]
%
% AUTHOR:
%   Turcu Alex Cristian, 2023, MATLAB, plotOrbit.m
% ------------------------------------------------------------------------

%% Arguments set

default_deltaT = -1;
default_n = 10000;
default_mu = astroConstants(13);
default_color = 'b';

checkVectorSize = @(x) isnumeric(x) && isvector(x) && isequal(length(x),3);
checkPositiveScalar = @(x) isnumeric(x) && isscalar(x) && (x>0);
checkPositiveInteger = @(x) isnumeric(x) && isscalar(x) && isinteger(x) && (x>0);

p = inputParser;
addRequired(p, 'r0', checkVectorSize)
addRequired(p, 'v0', checkVectorSize)
addOptional(p, 'deltaT', default_deltaT, checkPositiveScalar)
addOptional(p, 'n', default_n, checkPositiveInteger)
addParameter(p, 'mu', default_mu, checkPositiveScalar)
addParameter(p, 'color', default_color, @ischar)

parse(p, r0, v0, varargin{:})

pr = p.Results;
deltaT = pr.deltaT;
n = pr.n;
mu = pr.mu;
color = pr.color;

%% Time span set

if deltaT < 0
    N = -deltaT;
    a = 1/( 2/norm(r0) - dot(v0,v0)/mu );
    T = 2*pi*sqrt( a^3/mu );
    deltaT = T * N;
end

tspan = linspace( 0, deltaT, n );

%% State vector assembly

if size(r0, 2) ~= 1
    r0 = r0.';
end

if size(v0, 2) ~= 1
    v0 = v0.';
end

y0 = [ r0; v0 ];

%% ODE solver

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp(t, y, mu), tspan, y0, options );

%% Plot

plot3( Y(:,1), Y(:,2), Y(:,3), 'Color', color, 'LineWidth',2 )

%% Output

if nargout ~= 0
    t = T;
    r = Y(:, 1:3);
    v = Y(:, 4:6);
end