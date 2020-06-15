%% Autocorrelation analysis
% load the intensity file (optional if data is already in the workspace)
% and compute the autocorrelation (corr, cov, corr2, g, G)

load_status = 1;
ns = 60;
if load_status == 1
fileID = fopen('intensity_T10_ns100.bin');
intensity_vec = fread(fileID, ns*ntime, 'double');
fclose(fileID);

Intensity_dat = reshape(intensity_vec, [ntime, ns]);
clear intensity_vec
end
tic
%% 
ntau = 200;
t_dat = (1:ntime)*dt;
tau_dat = t_dat(1:ntau);

% nt - the number of time points used to for temporal average
nt = ntime - ntau;
%--------------------------------
% Override nt
%nt = 10^4;
%---------------------------------

%%
Corr_dat = zeros(ntau, ns);
Corr3_dat = zeros(ntau, ns);
G_dat = zeros(ntau, ns);
g_dat = zeros(ntau, ns);
I_avg = zeros(1, ns);

tic
for j = 1:ns
    I_avg(j) = mean(Intensity_dat(1:(nt+ntau),j));
    Intensity_1 = Intensity_dat(1:nt, j);
    meanI_1 = mean(Intensity_1);
    for k = 1:ntau
        Intensity_k = Intensity_dat(k:(k+nt-1),j);
        meanI_k = mean(Intensity_k);
        Corr_dat(k,j) = Intensity_1'*Intensity_k/nt;
        Corr3_dat(k,j) =  Intensity_1.^2'*Intensity_k/nt;
        G_dat(k,j) = Corr_dat(k,j)/(meanI_1*meanI_k) - 1;
        g_dat(k,j) = G_dat(k,j)*(meanI_1*meanI_k)/(std(Intensity_1)*std(Intensity_k));
    end
end
toc

clear Intensity_dat

Ne_ig = 1./G_dat(1,:);
F0_ig = I_avg*sqrt(8)./Ne_ig;

Corr_dat = Corr_dat';
Corr3_dat = Corr3_dat';
G_dat = G_dat';
g_dat = g_dat';