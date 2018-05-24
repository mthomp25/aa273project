% AA 273, Spring 2018
%
% 5/23/18
%
% Final Project
%
clear variables
close all

addpath('functions');

%% Get truth state
dur = 86400*3;% [s] Run for 3 days
dt = 5; % [s]

x_true = Truth_sim(dur, dt);