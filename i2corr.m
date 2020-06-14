%% Autocorrelation analysis
% load the intensity file (optional if data is already in the workspace)
% and compute the autocorrelation (corr, cov, corr2, g, G)

load_status = 1;
if load_status == 1
fileID = fopen('intensity.bin');
intensity_vec = fread(fileID, ns*ntime, 'double');
fclose(fileID);

Intensity_dat = reshape(intensity_vec, [ntime, ns])';
clear intensity_vec
end
tic
%% manually specify the number of tau points
ntau = 200;

tau_dat = t_dat(1:ntau);

% nt - the number of time points used to estimate autocorrlation E[I(t)I(t+tau)]
nt = ntime - ntau;

%--------------------------------
% Override nt
nt = 10^4;
%ns = 5;
%Intensity_dat = Intensity_dat(1:ns, 1:ntime);
%---------------------------------

%%
G_dat = zeros(ns, ntau);
Corr_dat = zeros(ns, ntau);
Corr3_dat = zeros(ns, ntau);
%Cov_dat = zeros(ns, ntau);
g2_dat = zeros(ns, ntau);
%Var_dat = zeros(ns, 1);
%I_dat = zeros(ns,1);

% Truncate each intensity siganl to the length of nt
Intensity_1 = Intensity_dat(1:ns, 1:nt); % ns by nt matrix
mean_Intensity_1 = mean(Intensity_1, 2); % ns by 1 column vector
var_1 = var(Intensity_1, 0, 2); % ns by 1 column vector
for k = 1:ntau
    Intensity_k = Intensity_dat(1:ns, k:(k+nt-1));
    mean_Intensity_k = mean(Intensity_k, 2);  
    var_k = var(Intensity_1, 0, 2); % ns by 1 column vector
    % dot(A, B, 2) computes the inner product of the rows of matrix A, B
    % where A, B are m by n matrices, dot(A, B, 2) is a m by 1 vector
    Corr_dat(:,k) = dot(Intensity_1, Intensity_k, 2)/nt;
    %Cov_dat(:,k) = Corr_dat(:,k) - mean_Intensity_1 .* mean_Intensity_k;
    Corr3_dat(:,k) = dot(Intensity_1.^2, Intensity_k, 2)/nt;
    G_dat(:,k) = Corr_dat(:,k) ./ (mean_Intensity_1 .* mean_Intensity_k) - 1;
    g2_dat(:,k) = Cov_dat(:,k) ./ sqrt(var_1 .* var_k);
end
I_dat = mean(Intensity_dat(1:ns, 1:nt),2);
Var_dat = var(Intensity_dat(1:ns,1:nt),0, 2);
% initial guess
% G(0) = 1/Ne, E[I] = F0*Ne/sqrt(8)
Ne_ig = 1./G_dat(:,1); % ns by 1 vector
F0_ig = I_dat*sqrt(8)./Ne_ig; % ns by 1 vector

toc