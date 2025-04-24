run("../Parameter_setting.m");
N_num = 5e2;
sx_vals = linspace(-Lx/2, Lx/2, N_num);
sy_vals = linspace(-Ly/2, Ly/2, N_num);

dlt_list = linspace(0,1,101);

cd ..;

for t = 1:length(dlt_list)
    t
    dlt = dlt_list(t);
    [~,~,j_now] = K_delta(dlt,1,N_num);
    save(['Current_density/Current_data/current_dlt=',num2str(dlt,'%.3f'),'.mat'],'j_now','-mat');
end

%contourf(sx_vals,sy_vals,magnitudes,60,'LineColor','none');

cd Current_density;
