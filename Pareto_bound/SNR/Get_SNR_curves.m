run("../../Parameter_setting.m");
load("../../Kmatrix_values.mat");

SNR_dB = -10:0.5:40;
SNR = 10.^(SNR_dB/10);

% sensinhg-opt

[eigenvec_r,eigenval_r] = eig(Krr);
[~,idx] = max(diag(eigenval_r));
k_r = eigenvec_r(:,idx);
mu_r = eigenval_r(idx,idx);

lambda_s = 6*RCS*gain*SNR*mu_r;
PDs = marcumq(sqrt(real(lambda_s)),7.4621,9);
Rs = real(log2(1+k_r'*Kur'*Kur*k_r/mu_r*SNR));

% communication-opt

[eigenvec_u,eigenval_u] = eig(Kuu);
[~,idx] = max(diag(eigenval_u));
k_u = eigenvec_u(:,idx);
mu_u = eigenval_u(idx,idx);

lambda_c = 6*RCS*gain*k_u'*Kur*Kur'*k_u/mu_u*SNR;
PDc = marcumq(sqrt(real(lambda_c)),7.4621,9);
Rc = log2(1+mu_u*SNR);

save('SNR_data.mat','-mat');