function [eu_opt,er_opt,j_delta] = K_delta(Krr,Kuu,Kur,delta,current_flag)
%run("Parameter_setting.m");
Pt = 0.1;
K_delta_mat = [delta*Kuu,sqrt(delta*(1-delta))*Kur;
    sqrt(delta*(1-delta))*Kur',(1-delta)*Krr];

[eigenvectors, eigenvalues] = eig(K_delta_mat);
[~,idx] = max(diag(eigenvalues));
k_delta = eigenvectors(:,idx);

mu = eigenvalues(idx,idx);
%%% optimal e fields
er_opt = sqrt(Pt/mu)*([sqrt(delta)*Kur',sqrt(1-delta)*Krr])*k_delta;
eu_opt = sqrt(Pt/mu)*([sqrt(delta)*Kuu,sqrt(1-delta)*Kur])*k_delta;

%% solve for optimal current
if current_flag == 1
    N_num = 1e2;
    sx_vals = linspace(-Lx/2, Lx/2, N_num);
    sy_vals = linspace(-Ly/2, Ly/2, N_num);
    
    j_delta = cell(N_num,N_num);
    
    % for i = 1:N_num
    %     for j = 1:N_num
    %         s = [sx_vals(i);sy_vals(j);0];
    %         j_delta{i,j} = 1j*k*Z0*[sqrt(delta)*Green_function(k,u-s)',sqrt(1-delta)*Green_function(k,r-s)']*k_delta;       %%% jKZ0 is just a constant scaling
    %     end
    % end
    for j = 1:N_num % x-coordinate
        for i = 1:N_num % y-coordinate
            s = [sx_vals(j);sy_vals(i);0];
            j_delta{i,j} = sqrt(Pt)/sqrt(mu)*1j*k*Z0*[sqrt(delta)*Green_function(k,u-s)',sqrt(1-delta)*Green_function(k,r-s)']*k_delta;
        end
    end
else
    j_delta = 0;
end

end