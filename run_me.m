% Run this script to run the numerical experiments and recreate the
% figures. It takes roughly 2 minutes. Run plotter_run_me.m to only reproduce the plots. 

clc
clear

rng(0)
addpath('Other')
addpath('results')

tic

pretty_function_example() % Figure 1
matern_example() % Figure 2
gaussian_process_example() % Figure 3
inverse_example() % Figure 4

toc