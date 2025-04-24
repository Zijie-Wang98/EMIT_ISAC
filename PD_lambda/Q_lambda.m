close all
color_blue = [0 61 124;
    98 159 208;
    146 195 228;
    194 230 247;
]/255;
LW = 1.5;
mark = {'-',':','--','-.'}; 


b = [7.0135 7.244 7.4621 7.6698];
b = flip(b);
PFA = 10.^[-5.5 -5 -4.5 -4];


lambda = 1:100;

figure;

for i = 1:length(b)
    P_D = marcumq(sqrt(lambda),b(i),9);
    semilogx(lambda,P_D,'LineStyle',mark{i},'LineWidth',LW,'Color',color_blue(i,:),'DisplayName',['$P_{\rm FA}=10^{',num2str(log10(PFA(i)),'%.1f'),'}$']);
    hold on;
end

grid on;
box on;

xlabel('$\lambda$','Interpreter','latex');
ylabel('$P_{\rm D}$','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',7);
legend('Interpreter','latex','Location','northwest','FontSize',7.5);

set(gcf,"Units",'centimeters',"Position",[5,5,6,4],'PaperSize',[6,4]);
set(gcf,"Renderer",'painters');

print(gcf,'../../Paper_writing/fig/PD_lambda','-dpdf');