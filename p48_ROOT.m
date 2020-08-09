% Michael R. Buche
% MAE 5730 - Problem 48
% Final Computation Project
% Root file

clear; close all; clc;

% Method of obtaing equations of motion
p.method = 'minimal_AMB';
%p.method = 'minimal_Lagrange';
%p.method = 'maximal';

% Number of links
p.N = 16;
p.l = linspace(1,2,p.N);
% Timespan
tspan = [0 6];

% Options for simple, automatically set parameters
p.simple_system = true;
p.nice_initial_conditions = false;
p.torques = false;
p.springs = false;
p.friction = false;

% Parameters setup
if p.simple_system
    p = p48_simple_system(p);
else
    % define parameters yourself here!
end

% Plotting fontsize
p.fontsize = 16;

% ODE45 accuracy settings
opts.reltol = 1e-3;
opts.abstol = 1e-3;

%% Solve for the motion
[t,x_G,y_G,x_end,y_end,theta,theta_dot,v_G] = p48_solve(tspan,p);

%% Option to load saved runs and avoid waiting for solver
cd 'Saved Workspaces'
save 16-pendulum.mat  ;
%load('minimal_AMB-chaotic-15_links-20_seconds.mat');
%load('minimal_AMB-chaotic-8_links-20_seconds_torques_springs_friction.mat');
%load('minimal_AMB-nice-6_links-20_seconds.mat');
%load('minimal_AMB-chaotic-6_links-20_seconds.mat');
%load('minimal_AMB-nice-5_links-20_seconds_point.mat');
%load('minimal_AMB-nice-3_links-20_seconds.mat');
%load('minimal_AMB-chaotic-3_links-20_seconds.mat');
%load('minimal_AMB-chaotic-3_links-40_seconds_springs_friction.mat')
%load('minimal_AMB-chaotic-3_links-20_seconds_springs.mat')
%load('minimal_Lagrange-nice-5_links-20_seconds_point.mat');
%load('minimal_Lagrange-nice-3_links-20_seconds.mat');
%load('minimal_Lagrange-chaotic-3_links-20_seconds.mat');
%load('minimal_Lagrange-chaotic-3_links-20_seconds_springs.mat')
%load('maximal-nice-3_links-20_seconds.mat');
%load('maximal-chaotic-3_links-20_seconds.mat');
cd ..

%% Plot the motion
p48_plot(t,x_end,y_end,theta,p)
    
%% Plot the changes in energy
p48_energy(t,x_G,v_G,theta,theta_dot,p);

%% Animate the motion
p.animation_timescale = 1;
p.timer = false; % F.Y.I. using this can make animation choppy
p.hinges = true;
pos = [900 550 560 420];
p48_animate(t,x_end,y_end,p,pos);

%% Create and save a gif
return % comment to allow gifs to be made
p.num_points = 25*tspan(2);
p.timer = true;
p.gif_file_name = 'name.gif';
p48_gif(t,x_end,y_end,p);