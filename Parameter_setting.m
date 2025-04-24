%%% physical constant
c = 3e8;
fc = 12e9;
lambda = c/fc;
k = 2*pi/lambda;
noise_var = 5.6e-3;

Z0 = 376.7;

%%% antenna constant
Lx = 30*15e-3;
Ly = 30*15e-3;
N = 30*30;

%%% system setup

Pt = 100e-3;
RCS = 0.1;   %%% RCS is sigma^2

u = [5;-5;15];
%u = [5;-20;15];
r = [10;5;30];
