%% Process the fitted D using different methods and T
%% read binary files D_method_T.bin and plot the bias and L2 error


% bias_table, error table = 
%-----------------------------------------
%      | g2  | G | Corr | Corr3 | Corr2n3
% ----------------------------------------
% T_i, i = 1, 2, ..., nT

bias_table = zeros(nT, 5);
error_table = zeros(nT, 5);
bias_errorbar = zeros(nT, 5);
error_lbar = zeros(nT, 5);
error_ubar = zeros(nT, 5);


D_corr_filenames = ["D_corr_T1.bin", "D_corr_Tp5.bin", "D_corr_Tp2.bin",...
    "D_corr_Tp1.bin", "D_corr_Tp05.bin", "D_corr_Tp02.bin", "D_corr_Tp01.bin"];
D_corr3_filenames = ["D_corr3_T1.bin", "D_corr3_Tp5.bin", "D_corr3_Tp2.bin",...
    "D_corr3_Tp1.bin", "D_corr3_Tp05.bin", "D_corr3_Tp02.bin", "D_corr3_Tp01.bin"];
D_G_filenames = ["D_G_T1.bin", "D_G_Tp5.bin", "D_G_Tp2.bin",...
    "D_G_Tp1.bin", "D_G_Tp05.bin", "D_G_Tp02.bin", "D_G_Tp01.bin"];
D_g2_filenames = ["D_g2_T1.bin", "D_g2_Tp5.bin", "D_g2_Tp2.bin",...
    "D_g2_Tp1.bin", "D_g2_Tp05.bin", "D_g2_Tp02.bin", "D_g2_Tp01.bin"];
D_2n3_filenames = ["D_2n3_T1.bin", "D_2n3_Tp5.bin", "D_2n3_Tp2.bin",...
    "D_2n3_Tp1.bin", "D_2n3_Tp05.bin", "D_2n3_Tp02.bin", "D_2n3_Tp01.bin"];

% g2
for i = 1:nT
    fileID = fopen(D_g2_filenames(i));
    D_fit = fread(fileID, 'double');
    fclose(fileID);
    
    diff = D_fit - D;
    diff2 = diff.^2;
    bias_table(i,1) = mean(diff);
    bias_errorbar(i,1) = std(diff)/sqrt(ns);
    error_table(i,1) = sqrt(mean(diff2));
    std_error2_bar = std(diff2)/sqrt(ns);
    error_ubar(i,1) = sqrt(mean(diff2)+std_error2_bar) - sqrt(mean(diff2));
    error_lbar(i,1) = sqrt(mean(diff2)) - sqrt(mean(diff2)-std_error2_bar);  
end

% G
for i = 1:nT
    fileID = fopen(D_G_filenames(i));
    D_fit = fread(fileID, 'double');
    fclose(fileID);
    
    diff = D_fit - D;
    diff2 = diff.^2;
    bias_table(i,2) = mean(diff);
    bias_errorbar(i,2) = std(diff)/sqrt(ns);
    error_table(i,2) = sqrt(mean(diff2));
    std_error2_bar = std(diff2)/sqrt(ns);
    error_ubar(i,2) = sqrt(mean(diff2)+std_error2_bar) - sqrt(mean(diff2));
    error_lbar(i,2) = sqrt(mean(diff2)) - sqrt(mean(diff2)-std_error2_bar);  
end

% Corr
for i = 1:nT
    fileID = fopen(D_corr_filenames(i));
    D_fit = fread(fileID, 'double');
    fclose(fileID);
    
    diff = D_fit - D;
    diff2 = diff.^2;
    bias_table(i,3) = mean(diff);
    bias_errorbar(i,3) = std(diff)/sqrt(ns);
    error_table(i,3) = sqrt(mean(diff2));
    std_error2_bar = std(diff2)/sqrt(ns);
    error_ubar(i,3) = sqrt(mean(diff2)+std_error2_bar) - sqrt(mean(diff2));
    error_lbar(i,3) = sqrt(mean(diff2)) - sqrt(mean(diff2)-std_error2_bar);  
end

% Corr3
for i = 1:nT
    fileID = fopen(D_corr3_filenames(i));
    D_fit = fread(fileID, 'double');
    fclose(fileID);
    
    diff = D_fit - D;
    diff2 = diff.^2;
    bias_table(i,4) = mean(diff);
    bias_errorbar(i,4) = std(diff)/sqrt(ns);
    error_table(i,4) = sqrt(mean(diff2));
    std_error2_bar = std(diff2)/sqrt(ns);
    error_ubar(i,4) = sqrt(mean(diff2)+std_error2_bar) - sqrt(mean(diff2));
    error_lbar(i,4) = sqrt(mean(diff2)) - sqrt(mean(diff2)-std_error2_bar);  
end

% Corr 2n3
for i = 1:nT
    fileID = fopen(D_2n3_filenames(i));
    D_fit = fread(fileID, 'double');
    fclose(fileID);
    
    diff = D_fit - D;
    diff2 = diff.^2;
    bias_table(i,5) = mean(diff);
    bias_errorbar(i,5) = std(diff)/sqrt(ns);
    error_table(i,5) = sqrt(mean(diff2));
    std_error2_bar = std(diff2)/sqrt(ns);
    error_ubar(i,5) = sqrt(mean(diff2)+std_error2_bar) - sqrt(mean(diff2));
    error_lbar(i,5) = sqrt(mean(diff2)) - sqrt(mean(diff2)-std_error2_bar);  
end

%%
rel_bias_table = bias_table/D*100;
rel_bias_errorbar = bias_errorbar/D*100;
rel_error_table = error_table/D*100;
rel_error_lbar = error_lbar/D*100;
rel_error_ubar = error_ubar/D*100;

figure
errorbar(T_dat, rel_error_table(:,1), rel_error_lbar(:,1), rel_error_ubar(:,1), '-o')
hold on
errorbar(T_dat, rel_error_table(:,2), rel_error_lbar(:,2), rel_error_ubar(:,2), '-*')
errorbar(T_dat, rel_error_table(:,3), rel_error_lbar(:,3), rel_error_ubar(:,3), '-x')
errorbar(T_dat, rel_error_table(:,4), rel_error_lbar(:,4), rel_error_ubar(:,4), '-+')
errorbar(T_dat, rel_error_table(:,5), rel_error_lbar(:,5), rel_error_ubar(:,5), '-s')
title('relative L2 error, ns = 200')
xlabel('T')
ylabel('relative L2 error, percent')
legend('g', 'G', 'Corr', 'Corr3', 'Corr2n3')
set(gca,'xscale','log', 'yscale', 'log')
grid on
saveas(gcf, 'rel_err_loglog_errorbar_2.png')

figure
errorbar(T_dat, rel_bias_table(:,1), rel_bias_errorbar(:,1), '-o')
hold on
errorbar(T_dat, rel_bias_table(:,2), rel_bias_errorbar(:,2), '-*')
errorbar(T_dat, rel_bias_table(:,3), rel_bias_errorbar(:,3),'-x')
errorbar(T_dat, rel_bias_table(:,4), rel_bias_errorbar(:,4),'-+')
errorbar(T_dat, rel_bias_table(:,5), rel_bias_errorbar(:,5), '-s')
title('relative bias, ns = 200')
xlabel('T')
ylabel('relative bias, percent')
legend('g', 'G', 'Corr', 'Corr3', 'Corr2n3')
set(gca,'xscale','log', 'yscale', 'log')
grid on
saveas(gcf, 'rel_bias_loglog_errorbar_2.png')

