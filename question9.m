%MCEN90018 - Advanced Fluid Dynamics - Assignment 2 - Q9
% This script answers Question 9 in the Assignment 2 sheet.

% Author: Francesco Mascadri
% Contact: fmascadri@student.unimelb.edu.au
% April 2022

%% Clear workspace
clear all;
clc;
close all;

%% Provided
g = 9.81;
w = 200000; %kg
Uinf_cruise = 250; %m/s
h_cruise = 10000; %m
rho_cruise = 0.4135; %kg/m3
AOA_cruise = 3; %degrees

rho_land = 1.225; %kg/m3
AOA_land = 10; %degrees

%% Found
CL_3 = 0.25; %est from Anderson graph (pg. 353)
CL_10 = 1.1;

%% Solve for wing area
S = (2*w*g)/(rho_cruise*Uinf_cruise^2*CL_3);

%% Solve for minimum landing speed
minLandingSpeed = sqrt((2*w*g)/(rho_land*S*CL_10));