% plot_actual_vs_est_cond.m
% Plot actual conductances versus the estimated
figure();
hold on;
subplot(2,1,1)
hold on
plot(that,gEhat,'-','Color',[0.4 0.4 1],'LineWidth',2);
plot(t,gE,'-k','LineWidth',2);
xlabel('time (ms)','FontSize',16);
ylabel(' g_E(t) (mS/cm^2)','FontSize',16);
legend('estimated','actual');
set(gca,'FontSize',14);
hold off;
 
subplot(2,1,2)
hold on
plot(that,gIhat,'-','Color',[1,0.4,0.6],'LineWidth',2);
plot(t,gI,'-k','LineWidth',2);
xlabel('time (ms)','FontSize',16);
ylabel('g_I(t) (mS/cm^2)','FontSize',16);
legend('estimated','actual');
set(gca,'FontSize',14);
hold off;

