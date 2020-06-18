%% Process the fitted D using different methods and T
%% read binary files D_method_T.bin and plot the bias and L2 error


% bias_table, error table = 
%-----------------------------------------
%      | g  | G | Corr | Corr3 | Corr2n3
% ----------------------------------------
% T_i, i = 1, 2, ..., nT

bias_table = zeros(nT, 5);
error_table = zeros(nT, 5);

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
