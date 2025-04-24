clc,clear all
run("Parameter_setting.m");
Kuu = Compute_K_integral(u,u);
Kur = Compute_K_integral(u,r);
Krr = Compute_K_integral(r,r);

sx = -Lx/2+(0.5+0:sqrt(N-1))*Lx/sqrt(N);
sy = -Ly/2+(0.5+0:sqrt(N-1))*Ly/sqrt(N);

%%% calculate \bar{g}(r)
gain = 0;
for i = 1:sqrt(N)
    for j = 1:sqrt(N)
        gain = gain+(4*pi*vecnorm([sx(i);sy(j);0]-r))^(-2);
    end
end

pho = max(eig(Kur))/sqrt(max(eig(Krr))*max(eig(Kuu)));

save('Kmatrix_values.mat',"Kuu","Krr","Kur","gain","pho",'-mat');

function I = Compute_K_integral(r1,r2)    
    run("Parameter_setting.m");
    N_num = 1e4;
    sx_vals = linspace(-Lx/2, Lx/2, N_num);
    sy_vals = linspace(-Ly/2, Ly/2, N_num);

    % calculate integration steps
    dx = Lx / (N_num - 1);
    dy = Ly / (N_num - 1);

    I = zeros(3, 3); % 3x3 matrix

    % calculate integrals
    for i = 1:N_num
        for j = 1:N_num
            s = [sx_vals(i); sy_vals(j); 0]; 
            G1 = Green_function(k,r1-s);
            G2 = Green_function(k,r2-s);
            integrand = G1 * G2'; 
            I = I + integrand * dx * dy;
        end
    end
    I = (k*Z0)^2*I;
end