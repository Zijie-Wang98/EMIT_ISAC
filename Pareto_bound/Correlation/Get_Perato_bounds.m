clc,clear all
run("../../Parameter_setting.m");

delta_list = 0:0.0001:1;
N_num = 3;
rho_list = zeros(N_num,1);
mags_eu = zeros(N_num,length(delta_list));
mags_er = mags_eu;
R = mags_eu;
PD = mags_eu;

for n = 1:N_num
    load(['Kmatrix_values',num2str(n),'.mat']);
    rho_list(n) = pho;
    for d = 1:length(delta_list)
        [eu,er,~] = K_delta(Krr,Kuu,Kur,delta_list(d),0);
        mags_eu(n,d) = vecnorm(eu)^2;
        mags_er(n,d) = vecnorm(er)^2;
    end
    gain
    R(n,:) = log2(1+mags_eu(n,:)/noise_var);
    lambda = (6*RCS*gain/noise_var)*mags_er(n,:);
    PD(n,:) = marcumq(sqrt(lambda),7.4621,9);
end

save('Pareto.mat','-mat');