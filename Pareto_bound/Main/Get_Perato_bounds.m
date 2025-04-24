clc,clear all
run("../../Parameter_setting.m");
load("../../Kmatrix_values.mat");
delta_list = 0:0.0001:1;
%delta_list = [linspace(0.1,0.1276,10) linspace(0.1276,0.1279,100) linspace(0.1279,5,5)];
mags_eu = zeros(1,length(delta_list));
mags_er = mags_eu;

cd ../..
for d = 1:length(delta_list)
    [eu,er,~] = K_delta(delta_list(d),0,0);
    mags_eu(d) = vecnorm(eu)^2;
    mags_er(d) = vecnorm(er)^2;
end
cd Pareto_bound/Main
R = log2(1+mags_eu/noise_var);


threshold0 = 7.6698; % Pfa = 1e-5.5
threshold1 = 7.4621; % Pfa = 1e-5
threshold2 = 7.244; % Pfa = 1e-4.5
threshold3 = 7.0135; % Pfa = 1e-4
Pfa_list = [10^(-5.5), 1e-5, 10^(-4.5), 1e-4];
lambda = (6*RCS*gain/noise_var)*mags_er;
PD0 = marcumq(sqrt(lambda),threshold0,9);
PD1 = marcumq(sqrt(lambda),threshold1,9);
PD2 = marcumq(sqrt(lambda),threshold2,9);
PD3 = marcumq(sqrt(lambda),threshold3,9);
