figure
h(1) = semilogy(E_SI_set_norm_all, 'LineWidth', 2);
hold on;
data_length = size(E_SI_set_norm_all, 1);
hold on;
h(2) = semilogy(E_DDOILC_norm, 'LineWidth', 2);
hold on;
h(3) = semilogy(E_SPSA2_all, 'LineWidth', 2);
hold on;
h(4) = semilogy(E_SPSA1_all, '--', 'LineWidth', 2, 'Color', 'green');
hold on;
h(5) = semilogy(E_SI_set_norm_all_reduction, 'LineWidth', 2, 'Color', 'black');
legend({'TRILC without trial reduction', 'DDOILC', 'SPSA2', 'SPSA1', 'TRILC with trial reduction'}, 'Interpreter','latex', 'location', 'southwest');
set(gca,'Children',[h(4) h(3) h(2) h(5) h(1)])
xlabel('Total trial number', 'Interpreter', 'latex'); 
ylabel('Tracking error $\frac{1}{2}\|\mathbf{E}\|^{2}$', 'Interpreter', 'latex');
set(get(gca,'XLabel'),'FontSize',15);
set(get(gca,'YLabel'),'FontSize',15); 
set(gcf,'unit','centimeters','position',[1 10 18 7]); 
grid minor