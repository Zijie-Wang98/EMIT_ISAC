clc,clear all
run("../../Parameter_setting.m");
r = [10;5;30];  %% case 1, rho = 0.0215;
r = [7;-10;30]; %% case 2, rho = 0.1994;
r = u;          %% case 3, rho = 1;
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
abs(pho)

save('Kmatrix_values3.mat',"Kuu","Krr","Kur","gain","pho",'-mat');

function I = Compute_K_integral(r1,r2)    
    % 生成积分网格
    run("../../Parameter_setting.m");
    N_num = 1e3;
    sx_vals = linspace(-Lx/2, Lx/2, N_num);
    sy_vals = linspace(-Ly/2, Ly/2, N_num);

    % 计算步长
    dx = Lx / (N_num - 1);
    dy = Ly / (N_num - 1);

    % 初始化积分结果
    I = zeros(3, 3); % G/G^H 是 3x3 矩阵
    cd ../../
    % 计算积分
    for i = 1:N_num
        for j = 1:N_num
            s = [sx_vals(i); sy_vals(j); 0]; % 取当前点 s
            G1 = Green_function(k,r1-s);
            G2 = Green_function(k,r2-s);
            integrand = G1 * G2'; % 计算 G*G^H
            I = I + integrand * dx * dy; % 乘以面积元素 dx*dy
        end
    end
    cd Pareto_bound/Correlation
    I = (k*Z0)^2*I;
end