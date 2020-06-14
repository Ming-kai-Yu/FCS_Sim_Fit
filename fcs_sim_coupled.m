%% fcs simulation driver

ns = 100;
fprintf("-Number of simulations: %d\n", ns);

dt = 1e-06;
fprintf("Time step dt: %g\n", dt);

T = 1e-01;
fprintf("-T: %g\n", T);
taumax = 2e-04;

D = 1.3e-10;
fprintf("-Diffusivity: %g\n", D);

r = 2e-7;
fprintf("Optical system r: %g\n", r);
l = 1e-6;
fprintf("Optical system l: %g\n", l);

br = 8*r;
fprintf("Box length: %g\n", 2*br);
bl = 8*l;
fprintf("Box height: %g\n", 2*bl);
conc = 0.02;
fprintf("Concentration: %g\n", conc);



% derived simulation parameters 
vol_sol = 8*br*br*bl;
avogadro = 6.0221409e+23;
npart_d = conc*avogadro*(1e-06)*1000*vol_sol;
npart_i = round(npart_d);
fprintf("Number of particles %d\n", npart_i);

vol_det = pi * sqrt(pi) * r*r*l;
neff = conc*avogadro*(1e-06)*1000*vol_det;
fprintf("-Effective number of particle in illuminatd region: %g\n", neff);

ntime = round((T+taumax)/dt); 
fprintf("-The number of time points: %d\n", ntime);

t_dat = (1:ntime)*dt;

%% Run the simulation
a = 3*r;
c = 3*l;
m = 10;
fprintf('m = %d.\n', m);

tic

intensity_acc_mat = zeros(ntime, ns);
intensity_eff_mat = zeros(ntime, ns);

for j = 1:ns
    %[intensity_acc_mat(:,j), intensity_eff_mat(:,j)] = get_intensity_coupled_ellip...
    %(dt, ntime, D, br, bl, npart_i, r, l, a, c, m);
    [intensity_acc_mat(:,j), intensity_eff_mat(:,j)] = get_intensity_coupled_new...
    (dt, ntime, D, br, bl, npart_i, r, l, m);
end
toc


%% Compare intensity
intensity_diff = intensity_acc_mat - intensity_eff_mat; % nt by ns matrix
intensity_diff_mean = mean(intensity_diff, 2);
intensity_diff_std = std(intensity_diff, 0, 2);

intensity_diff_sup_t = max(intensity_diff); 
%% sup_{t in [0,T]} [Ia(t) - Ie(t)]
fprintf('--------Comparing coupled Ie with Ia-------------\n')
fprintf('The mean of e = sup_{t in [0,T]} [Ie(t) - Ia(t)] is %g.\n',...
    mean(intensity_diff_sup_t));
fprintf('The std of sup_{t in [0,T]} [Ia(t) - Ie(t)] is %g.\n',...
    std(intensity_diff_sup_t));
fprintf('----The confidence interval for e is [%g, %g].\n',...
    mean(intensity_diff_sup_t) - std(intensity_diff_sup_t)/sqrt(ns), ...
    mean(intensity_diff_sup_t) + std(intensity_diff_sup_t)/sqrt(ns));
fprintf('The L2 error, i.e. sqrt of the mean of e^2 is %g.\n',...
    sqrt(mean(intensity_diff_sup_t.^2)));
fprintf('The std of e^2 is %g.\n', std(intensity_diff_sup_t.^2));
fprintf('The upper limit of sqrt(e^2_bar) is %g.\n',...
    sqrt(mean(intensity_diff_sup_t.^2) + std(intensity_diff_sup_t.^2)/sqrt(ns)));
fprintf('The lower limit of sqrt(e^2_bar) is %g.\n',...
    sqrt(mean(intensity_diff_sup_t.^2) - std(intensity_diff_sup_t.^2)/sqrt(ns)));

%%
% temporal average
intensity_acc_avg = mean(intensity_acc_mat); % 1 by ns vector
intensity_eff_avg = mean(intensity_eff_mat); % 1 by ns vector
intensity_acc_avg = intensity_acc_avg'; % ns by 1 vector
intensity_eff_avg = intensity_eff_avg'; % ns by 1 vector

var_intensity_acc = var(intensity_acc_mat);
var_intensity_eff = var(intensity_eff_mat);
var_intensity_acc = var_intensity_acc';
var_intensity_eff = var_intensity_eff';

mean_intensity_acc = mean(intensity_acc_mat, 2);
std_intensity_acc = std(intensity_acc_mat, 0, 2);
mean_intensity_eff = mean(intensity_eff_mat, 2);
std_intensity_eff = std(intensity_eff_mat, 0, 2);
%%
%{
figure
errorbar(t_dat(1:2000:ntime), intensity_diff_mean(1:2000:ntime),...
    intensity_diff_std(1:2000:ntime)/sqrt(ns))
xlabel('t')
ylabel('Ia - Ie')
%saveas(gcf,'intensity_diff_errorbar_sparse.png')

figure
plot(t_dat(1:50), intensity_accurate(1:50), '-o')
hold on
plot(t_dat(1:50), intensity_efficient(1:50), '-+')
plot(t_dat(1:50), intensity_accurate(1:50)-intensity_efficient(1:50), '-*')
xlabel('t')
ylabel('intensity')
legend('accurate', 'efficient', 'difference')
%saveas(gcf, 'sample_acc_eff_portion.png')

figure
errorbar(t_dat(1:2000:ntime), mean_intensity_acc(1:2000:ntime), ...
    std_intensity_acc(1:2000:ntime)/sqrt(ns), '-o')
hold on 
errorbar(t_dat(1:2000:ntime), mean_intensity_eff(1:2000:ntime), ...
    std_intensity_eff(1:2000:ntime)/sqrt(ns), '-x')
errorbar(t_dat(1:2000:ntime), intensity_diff_mean(1:2000:ntime),...
    intensity_diff_std(1:2000:ntime)/sqrt(ns))
xlabel('t')
ylabel('intensity')
legend('accurate', 'efficient', 'difference')
saveas(gcf, 'errorbar_acc_eff.png')
%}

%% Compare Corr
ntau = 200;
tau_dat = (1:ntau)*dt;
nt = ntime-ntau;
intensity_acc_dat = intensity_acc_mat';
intensity_eff_dat = intensity_eff_mat';
Corr_acc = zeros(ns, ntau);
Corr_eff = zeros(ns, ntau);

Intensity_1 = intensity_acc_dat(1:ns, 1:nt); % ns by nt matrix
for k = 1:ntau
    Intensity_k = intensity_acc_dat(1:ns, k:(k+nt-1));
    % dot(A, B, 2) computes the inner product of the rows of matrix A, B
    % where A, B are m by n matrices, dot(A, B, 2) is a m by 1 vector
    Corr_acc(:,k) = dot(Intensity_1, Intensity_k, 2)/nt;
end

Intensity_1 = intensity_eff_dat(1:ns, 1:nt); % ns by nt matrix
for k = 1:ntau
    Intensity_k = intensity_eff_dat(1:ns, k:(k+nt-1));
    % dot(A, B, 2) computes the inner product of the rows of matrix A, B
    % where A, B are m by n matrices, dot(A, B, 2) is a m by 1 vector
    Corr_eff(:,k) = dot(Intensity_1, Intensity_k, 2)/nt;
end

%%
corr_diff = Corr_acc - Corr_eff;
corr_diff_mean = mean(corr_diff);
corr_diff_std = std(corr_diff);

corr_acc_mean = mean(Corr_acc);
corr_eff_mean = mean(Corr_eff);

%{
figure
errorbar(tau_dat, corr_diff_mean, corr_diff_std)
xlabel('t')
ylabel('Corr_a - Corr_e')
%saveas(gcf,'corr_errorbar.png')

figure
errorbar(tau_dat(1:4:ntau), corr_diff_mean(1:4:ntau), corr_diff_std(1:4:ntau))
xlabel('t')
ylabel('Corr_a - Corr_e')
%saveas(gcf,'corr_errorbar_sparse.png')

figure
errorbar(tau_dat(1:4:ntau), corr_diff_mean(1:4:ntau), corr_diff_std(1:4:ntau)/sqrt(ns))
xlabel('t')
ylabel('Corr_a - Corr_e')
saveas(gcf,'mean_corr_errorbar_sparse.png')
figure
plot(tau_dat, corr_acc_mean)
hold on
plot(tau_dat, corr_eff_mean)
plot(tau_dat, corr_diff_mean)

%%
corr_diff_sup_t = max(corr_diff, [], 2);
fprintf('--------Comparing coupled Corr_e with Corr_a -------------\n')
fprintf('The mean of e = sup_{t in [0,T]} [Corr_e(t) - Corr_a(t)] is %g.\n',...
    mean(corr_diff_sup_t));
fprintf('The std of e is %g.\n',...
    std(corr_diff_sup_t));
fprintf('The confidence interval for e is [%g, %g].\n',...
    mean(corr_diff_sup_t) - std(corr_diff_sup_t)/sqrt(ns), ...
    mean(corr_diff_sup_t) + std(corr_diff_sup_t)/sqrt(ns));
fprintf('The L2 error, i.e. sqrt of the mean of e^2 is %g.\n',...
    sqrt(mean(corr_diff_sup_t.^2)));
fprintf('The std of e^2 is %g.\n', std(corr_diff_sup_t.^2));
fprintf('The upper limit of sqrt(e^2_bar) is %g.\n',...
    sqrt(mean(corr_diff_sup_t.^2) + std(corr_diff_sup_t.^2)/sqrt(ns)));
fprintf('The lower limit of sqrt(e^2_bar) is %g.\n',...
    sqrt(mean(corr_diff_sup_t.^2) - std(corr_diff_sup_t.^2)/sqrt(ns)));
%}
%% Compare fitted D
% k scaling factor, in order to make D and Ne, F0 at the same scale
k = 10^10;

Corr_model = @(params, tau_dat) params(3)^2 * ( params(1) * ...
    (1/8./(1+4*params(2)*tau_dat/(k*r^2))) .* (1./sqrt(1+4*params(2)*tau_dat/(k*l^2))) ...
    + params(1)^2/8);
% where params(1) is Neff, parmas(2) is D, and params(3) is F0

%% Initial Guess
Ne_ig_acc = intensity_acc_avg.^2./var_intensity_acc;
F0_ig_acc = sqrt(8)*var_intensity_acc./intensity_acc_avg;
Ne_ig_eff = intensity_eff_avg.^2./var_intensity_eff;
F0_ig_eff = sqrt(8)*var_intensity_eff./intensity_eff_avg;
%%
D_true = 1.3e-10;
params_fit_Corr_acc = zeros(ns, 3);
params_fit_Corr_eff = zeros(ns, 3);


fprintf('---------------------------------\n');
fprintf('Fitting Corr accurate.\n');
for i = 1:ns
    Corr_obs = Corr_acc(i, :);
    % params0 = [Ne_ig, D_ig, F0_ig];
    params0 = [Ne_ig_acc(i), 1, F0_ig_acc(i)];
    % f1 maps params to ntau dimensional vector
    f1 = @(params)  Corr_obs - Corr_model(params, tau_dat);
    options = optimoptions('lsqnonlin',...
        'Display', 'off', 'TolX', 1e-6);
    params_fit_Corr_acc(i,:) = lsqnonlin(f1, params0, [0,0,0], [], options);
end
D_fit_Corr_acc = params_fit_Corr_acc(:,2)/k;

%%
fprintf('T = %g.\n', nt*dt);
%fprintf('True parameters: Neff %g, D %g.\n', params_true(1), params_true(2));
fprintf('The mean of the error for fitted D is %g.\n', mean(D_fit_Corr_acc)-D_true);
fprintf('The relative bias is %g percent.\n', (mean(D_fit_Corr_acc)-D_true)/D_true*100);
fprintf('The std of the fitted D is %g.\n', std(D_fit_Corr_acc));
fprintf('The confidence interval for the error is [%g, %g].\n',...
    mean(D_fit_Corr_acc) - D_true - std(D_fit_Corr_acc)/sqrt(ns), ...
    mean(D_fit_Corr_acc) - D_true + std(D_fit_Corr_acc)/sqrt(ns));
fprintf('The relative L2 error is %g percent.\n', ...
    sqrt((mean(D_fit_Corr_acc)-D_true)^2 + std(D_fit_Corr_acc)^2)/D_true*100);
fprintf('The total number of trajectories is %g.\n', ns);
%%
fprintf('---------------------------------\n');
fprintf('Fitting Corr efficient.\n');
for i = 1:ns
    Corr_obs = Corr_eff(i, :);
    % params0 = [Ne_ig, D_ig, F0_ig];
    params0 = [Ne_ig_eff(i), 1, F0_ig_eff(i)];
    % f1 maps params to ntau dimensional vector
    f1 = @(params)  Corr_obs - Corr_model(params, tau_dat);
    options = optimoptions('lsqnonlin',...
        'Display', 'off', 'TolX', 1e-6);
    params_fit_Corr_eff(i,:) = lsqnonlin(f1, params0, [0,0,0], [], options);
end
D_fit_Corr_eff = params_fit_Corr_eff(:,2)/k;

%%
fprintf('T = %g.\n', nt*dt);
%fprintf('True parameters: Neff %g, D %g.\n', params_true(1), params_true(2));
fprintf('The mean of the error for fitted D is %g.\n', mean(D_fit_Corr_eff)-D_true);
fprintf('The relative bias is %g percent.\n', (mean(D_fit_Corr_eff)-D_true)/D_true*100);
fprintf('The std of the fitted D is %g.\n', std(D_fit_Corr_eff));
fprintf('The relative L2 error is %g percent.\n', ...
    sqrt((mean(D_fit_Corr_eff)-D_true)^2 + std(D_fit_Corr_eff)^2)/D_true*100);
fprintf('The total number of trajectories is %g.\n', ns);
%%
D_diff = D_fit_Corr_eff - D_fit_Corr_acc;
mean_D_diff = mean(D_diff);
var_D_diff = var(D_diff);

fprintf('--------Comparing coupled D_hat_e with D_hat_a-------------\n')
fprintf('The mean of d = D_hat_e - D_hat_a is %g.\n',...
    mean(D_diff));
fprintf('The std of d is %g.\n',...
    std(D_diff));
fprintf('The confidence interval for e is [%g, %g].\n',...
    mean(D_diff) - std(D_diff)/sqrt(ns), ...
    mean(D_diff) + std(D_diff)/sqrt(ns));
fprintf('The L2 error, i.e. sqrt of the mean of d^2 is %g.\n',...
    sqrt(mean(D_diff.^2)));
fprintf('The std of d^2 is %g.\n', std(D_diff.^2));
upper_l2 = sqrt(mean(D_diff.^2) + std(D_diff.^2)/sqrt(ns));
lower_l2 = sqrt(mean(D_diff.^2) - std(D_diff.^2)/sqrt(ns));
fprintf('The upper limit of sqrt(d^2_bar) is %g.\n',...
    upper_l2);
fprintf('The lower limit of sqrt(d^2_bar) is %g.\n',...
    lower_l2);
fprintf('----The confidence interval for relative error is [%g, %g].\n',...
    lower_l2/D, upper_l2/D);
    