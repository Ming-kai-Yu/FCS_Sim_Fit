%% Comparison of bias and error of various fit

% bias_table, error table = 
%-----------------------------------------
%      | g  | G | Corr | Corr3 | Corr2n3
% ----------------------------------------
% T_i, i = 1, 2, ..., nT
init_status = 0;
if init_status == 1
bias_table = zeros(7, 5);
error_table = zeros(7, 5);
bias_errorbar = zeros(7, 5);
error_lbar = zeros(7, 5);
error_ubar = zeros(7, 5);
end
%%
i = 7;
T_dat(i) = 0.01; 
D_true = 1.3e-10;

error = g_G_Corr_Corr3_Corr2n3_N2Tp01 - D_true;
bias = mean(error);
e_std = std(error);
e_bar_std = e_std/sqrt(ns);
bias_lower = bias - e_bar_std;
bias_upper = bias + e_bar_std;

e2 = error.^2;
e2_bar = mean(e2);
l2 = sqrt(e2_bar);
e2_std = std(e2);
e2_bar_std = e2_std/sqrt(ns);
error_lower = sqrt(e2_bar - e2_bar_std);
error_upper = sqrt(e2_bar + e2_bar_std);

bias_table(i,:) = bias;
error_table(i,:) = l2;
bias_errorbar(i,:) = e_bar_std;
error_lbar(i,:) = l2 - error_lower;
error_ubar(i,:) = error_upper - l2;

%%
rel_bias_table = bias_table/D*100;
rel_error_table = error_table/D*100;
rel_error_lbar = error_lbar/D*100;
rel_error_ubar = error_ubar/D*100;
rel_bias_errorbar = bias_errorbar/D*100;


%%
figure
errorbar(T_dat, rel_error_table(:,1), rel_error_lbar(:,1), rel_error_ubar(:,1), '-o')
hold on
errorbar(T_dat, rel_error_table(:,2), rel_error_lbar(:,2), rel_error_ubar(:,2), '-*')
errorbar(T_dat, rel_error_table(:,3), rel_error_lbar(:,3), rel_error_ubar(:,3), '-x')
errorbar(T_dat, rel_error_table(:,4), rel_error_lbar(:,4), rel_error_ubar(:,4), '-+')
%errorbar(T_dat, rel_error_table(:,5), rel_error_lbar(:,5), rel_error_ubar(:,5))
xlabel('T')
ylabel('relative L2 error, percent')
legend('g', 'G', 'Corr', 'Corr3', 'Corr2n3')
set(gca,'xscale','log', 'yscale', 'log')
grid on
saveas(gcf, 'relerr_loglog_errorbar.png')



%%
figure
loglog(T_dat, rel_bias_table(:,1), '-o')
hold on
loglog(T_dat, rel_bias_table(:,2), '-*')
loglog(T_dat, rel_bias_table(:,3), '-x')
loglog(T_dat, rel_bias_table(:,4), '-+')
loglog(T_dat, rel_bias_table(:,5))
xlabel('T')
ylabel('relative bias, percent')
legend('g', 'G', 'Corr', 'Corr3', 'Corr2n3')
grid on
saveas(gcf, 'rel_bias_loglog_1.png')

figure
errorbar(T_dat, rel_bias_table(:,1), rel_bias_errorbar(:,1), '-o')
hold on
errorbar(T_dat, rel_bias_table(:,2), rel_bias_errorbar(:,2), '-*')
errorbar(T_dat, rel_bias_table(:,3), rel_bias_errorbar(:,3),'-x')
errorbar(T_dat, rel_bias_table(:,4), rel_bias_errorbar(:,4),'-+')
%errorbar(T_dat, rel_error_table(:,5), rel_error_lbar(:,5), rel_error_ubar(:,5))
xlabel('T')
ylabel('relative bias, percent')
legend('g', 'G', 'Corr', 'Corr3', 'Corr2n3')
set(gca,'xscale','log', 'yscale', 'log')
grid on
saveas(gcf, 'relbias_loglog_errorbar.png')



%%
figure
loglog(T_dat, rel_error_table(:,1), '-o')
hold on
loglog(T_dat, rel_error_table(:,2), '-*')
loglog(T_dat, rel_error_table(:,3), '-x')
loglog(T_dat, rel_error_table(:,4), '-+')
loglog(T_dat, rel_error_table(:,5))
xlabel('T')
ylabel('relative L2 error, percent')
legend('g', 'G', 'Corr', 'Corr3', 'Corr2n3')
grid on
saveas(gcf, 'rel_error_loglog_1.png')