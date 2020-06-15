% Perform different methods of fitting
% written by Mingkai Yu, 03/09/20
% edited by Mingkai Yu 06/10/20
tic
%% load data

% manually specify optical system parameters r and l
r = 2e-7;
l = 1e-6;

% manually input the ture parameters
Neff_true = 2.6826604619468237; % dimensionless
D_true = 1.3e-10; % meter^2/second
F0 = 1;
params_true = [Neff_true, D_true, F0];
tau_r = r^2/(4*D_true);


%% Models 
k = 10^10;
% Gfun_model is the scaled autocorrelation function 
% since D is on the 10^(-10) scale while Neff is on the scale of 1
Gfun_model = @(params, tau_dat) 1/params(1) ...
    * (1./(1+4*params(2)*tau_dat/(k*r^2))) .* (1./sqrt(1+4*params(2)*tau_dat/(k*l^2)));
% where params(1) is Neff, the effective number of particles in illuminated volume
% and params(2) is D, the diffusivity

Cov_model = @(params, tau_dat) params(1) ...
    * (1./(1+4*params(2)*tau_dat/(k*r^2))) .* (1./sqrt(1+4*params(2)*tau_dat/(k*l^2)));
% where params(1) is F0^2*Neff/8, params(2) is D

Corr_model = @(params, tau_dat) params(3)^2 * ( params(1) * ...
    (1/8./(1+4*params(2)*tau_dat/(k*r^2))) .* (1./sqrt(1+4*params(2)*tau_dat/(k*l^2))) ...
    + params(1)^2/8);
% where params(1) is Neff, parmas(2) is D, and params(3) is F0

Corr3_model = @(params, tau_dat) params(3)^3 * ( params(1)/sqrt(8) * ...
    1./(3+16*params(2)*tau_dat/(k*r^2)) .* 1./sqrt(3+16*params(2)*tau_dat/(k*l^2)) ...
    + params(1)^3/(8*sqrt(8)) + params(1)^2/(8*sqrt(8))...
    + params(1)^2/(8*sqrt(8))*2./(1+4*params(2)*tau_dat/(k*r^2)) ...
    .* (1./sqrt(1+4*params(2)*tau_dat/(k*l^2))) );
% where params(1) is Neff, parmas(2) is D, and params(3) is F0

Corr_2n3_model = @(params, tau_dat) ...
    Corr_model(params, tau_dat) + Corr3_model(params, tau_dat);

M2bar_model = @(param, tau_dat) ( param(1) * ...
    (1/8./(1+4*param(2)*tau_dat/(k*r^2))) .* (1./sqrt(1+4*param(2)*tau_dat/(k*l^2))) ...
    + param(1)^2/8);
M3bar_model = @(param, tau_dat) ( param(1)/sqrt(8) * ...
    1./(3+16*param(2)*tau_dat/(k*r^2)) .* 1./sqrt(3+16*param(2)*tau_dat/(k*l^2)) ...
    + param(1)^3/(8*sqrt(8)) + param(1)^2/(8*sqrt(8))...
    + param(1)^2/(8*sqrt(8))*2./(1+4*param(2)*tau_dat/(k*r^2)) ...
    .* (1./sqrt(1+4*param(2)*tau_dat/(k*l^2))) );

g2_model = @(params, tau_dat) (1./(1+4*params*tau_dat/(k*r^2))) ...
    .* (1./sqrt(1+4*params*tau_dat/(k*l^2)));



%% Fitting G
group_sz = 1;
num_group = ns;

G_status = 1;
if G_status == 1
    params_fit = zeros(ns, 2);
    fprintf('---------------------------------\n');
    fprintf('Fitting G.\n');
    for i = 1:ns
        G_obs = G_dat(i, :);
        params0 = [Ne_ig(i), 1];
        % f1 maps params to ntau dimensional vector
        f1 = @(params)  G_obs - Gfun_model(params, tau_dat);
        options = optimoptions('lsqnonlin',...
            'Display', 'off', 'TolX', 1e-6);
        params_fit(i,:) = lsqnonlin(f1, params0, [0, 0], [], options);
    end
D_fit_G = params_fit(:,2)/k;


fprintf('T = %g.\n', nt*dt);
fprintf('True parameters: Neff %g, D %g.\n', params_true(1), params_true(2));
fprintf('The mean of the error for fitted D is %g.\n', mean(D_fit_G)-D_true);
fprintf('The relative error is %g percent.\n', (mean(D_fit_G)-D_true)/D_true*100);
fprintf('The std of the fitted D is %g.\n', std(D_fit_G));
fprintf('The relative L2 error is %g percent.\n', ...
    sqrt((mean(D_fit_G)-D_true)^2 + std(D_fit_G)^2)/D_true*100);
fprintf('The total number of trajectories is %g.\n', ns);
fprintf('The number of least square fit performed is %d.\n', num_group);
end
%% fit Corr
corr_status = 1;
if corr_status == 1
    params_fit_Corr = zeros(ns, 3);
    fprintf('---------------------------------\n');
    fprintf('Fitting Corr.\n');
    for i = 1:ns
        Corr_obs = Corr_dat(i, :);
        params0 = [Ne_ig(i), 1, F0_ig(i)];
        % f1 maps params to ntau dimensional vector
        f1 = @(params)  Corr_obs - Corr_model(params, tau_dat);
        options = optimoptions('lsqnonlin',...
            'Display', 'off', 'TolX', 1e-6);
        params_fit_Corr(i,:) = lsqnonlin(f1, params0, [0,0,0], [], options);
    end
D_fit_Corr = params_fit_Corr(:,2)/k;


fprintf('T = %g.\n', nt*dt);
fprintf('True parameters: Neff %g, D %g.\n', params_true(1), params_true(2));
fprintf('The mean of the error for fitted D is %g.\n', mean(D_fit_Corr)-D_true);
fprintf('The relative error is %g percent.\n', (mean(D_fit_Corr)-D_true)/D_true*100);
fprintf('The std of the fitted D is %g.\n', std(D_fit_Corr));
fprintf('The relative L2 error is %g percent.\n', ...
    sqrt((mean(D_fit_Corr)-D_true)^2 + std(D_fit_Corr)^2)/D_true*100);
fprintf('The total number of trajectories is %g.\n', ns);
fprintf('The number of least square fit performed is %d.\n', num_group);
end
%% fit Corr3
corr3_status = 1;
if corr3_status == 1
    params_fit = zeros(ns, 3);
    fprintf('---------------------------------\n');
    fprintf('Fitting Corr3.\n');
    for i = 1:ns
        Corr3_obs = Corr3_dat(i, :);
        params0 = [Ne_ig(i), 1, F0_ig(i)];        
        % f1 maps params to ntau dimensional vector
        f1 = @(params)  Corr3_obs - Corr3_model(params, tau_dat);
        options = optimoptions('lsqnonlin',...
            'Display', 'off', 'TolX', 1e-6);
        params_fit(i,:) = lsqnonlin(f1, params0, [0,0,0], [], options);
    end
D_fit_Corr3 = params_fit(:,2)/k;


fprintf('T = %g.\n', nt*dt);
fprintf('True parameters: Neff %g, D %g.\n', params_true(1), params_true(2));
fprintf('The mean of the error for fitted D is %g.\n', mean(D_fit_Corr3)-D_true);
fprintf('The relative error is %g percent.\n', (mean(D_fit_Corr3)-D_true)/D_true*100);
fprintf('The std of the fitted D is %g.\n', std(D_fit_Corr3));
fprintf('The relative L2 error is %g percent.\n', ...
    sqrt((mean(D_fit_Corr3)-D_true)^2 + std(D_fit_Corr3)^2)/D_true*100);
fprintf('The total number of trajectories is %g.\n', ns);
fprintf('The number of least square fit performed is %d.\n', num_group);
end
%% Fitting g2
g2_status = 1;
if g2_status == 1
    fprintf('---------------------------------\n');
    fprintf('Fitting g2.\n');
    params_fit = zeros(num_group,1);
    params0 = 1;
    for i = 1:ns
        g2_obs = g2_dat(i,:);
        g2diff = @(params)  g2_obs - g2_model(params, tau_dat);
        options = optimoptions('lsqnonlin',...
            'Display', 'off', 'TolX', 1e-6);
        params_fit(i) = lsqnonlin(g2diff, params0, 0, [], options);
    end
D_fit_g2 = params_fit/k;


fprintf('T = %g.\n', nt*dt);
fprintf('True parameters: Neff %g, D %g.\n', params_true(1), params_true(2));
fprintf('The mean of the error for fitted D is %g.\n', mean(D_fit_g2)-D_true);
fprintf('The relative error is %g percent.\n', (mean(D_fit_g2)-D_true)/D_true*100);
fprintf('The std of the fitted D is %g.\n', std(D_fit_g2));
fprintf('The relative L2 error is %g percent.\n', ...
    sqrt((mean(D_fit_g2)-D_true)^2 + std(D_fit_g2)^2)/D_true*100);
fprintf('The total number of trajectories is %g.\n', ns);
fprintf('The number of least square fit performed is %d.\n', num_group);
end
%% Fitting Corr2 and Corr3
corr_2n3_status = 1;
if corr_2n3_status == 1
    params_fit = zeros(ns, 3);
    fprintf('---------------------------------\n');
    fprintf('Fitting Corr2 and Corr3.\n');
    for i = 1:ns
        Corr_obs = Corr_dat(i,:);
        Corr3_obs = Corr3_dat(i , :);
        % first step: fit Corr2 to get initial guess 
        params0 = params_fit_Corr(i,:);
        % second step: combined Corr2 and Corr3 fit
        f2 = @(params) [Corr_obs/params(3)^2 ...
        - M2bar_model([params(1), params(2)], tau_dat), ...
        Corr3_obs/params(3)^3 ...
        - M3bar_model([params(1), params(2)], tau_dat)];
        params_fit(i,:) = lsqnonlin(f2, params0, [0,0,0], [], options);
    end
    D_fit_Corr_2n3 = params_fit(:,2)/k;


fprintf('T = %g.\n', nt*dt);
fprintf('True parameters: Neff %g, D %g.\n', params_true(1), params_true(2));
fprintf('The mean of the error for fitted D is %g.\n', mean(D_fit_Corr_2n3)-D_true);
fprintf('The relative error is %g percent.\n', (mean(D_fit_Corr_2n3)-D_true)/D_true*100);
fprintf('The std of the fitted D is %g.\n', std(D_fit_Corr_2n3));
fprintf('The relative L2 error is %g percent.\n', ...
    sqrt((mean(D_fit_Corr_2n3)-D_true)^2 + std(D_fit_Corr_2n3)^2)/D_true*100);
fprintf('The total number of trajectories is %g.\n', ns);
fprintf('The number of least square fit performed is %d.\n', num_group);
end

%% Fitting Corr2 and Corr3, new normalization
corr_2n3_new_status = 0;
if corr_2n3_new_status == 1
    params_fit = zeros(ns, 3);
    fprintf('---------------------------------\n');
    fprintf('Fitting Corr2 and Corr3, normalize by (F*Ne)^k.\n');
    for i = 1:ns
        Corr_obs = Corr_dat(i,:);
        Corr3_obs = Corr3_dat(i , :);
        % first step: fit Corr2 to get initial guess 
        params0 = params_fit_Corr(i,:);
        % second step: combined Corr2 and Corr3 fit
        f2 = @(params) [Corr_obs/params(3)^2 ...
        - M2bar_model([params(1), params(2)], tau_dat), ...
        Corr3_obs/params(3)^3 ...
        - M3bar_model([params(1), params(2)], tau_dat)];
        params_fit(i,:) = lsqnonlin(f2, params0, [0,0,0], [], options);
    end
    D_fit_Corr_2n3 = params_fit(:,2)/k;


fprintf('T = %g.\n', nt*dt);
fprintf('True parameters: Neff %g, D %g.\n', params_true(1), params_true(2));
fprintf('The mean of the error for fitted D is %g.\n', mean(D_fit_Corr_2n3)-D_true);
fprintf('The relative error is %g percent.\n', (mean(D_fit_Corr_2n3)-D_true)/D_true*100);
fprintf('The std of the fitted D is %g.\n', std(D_fit_Corr_2n3));
fprintf('The relative L2 error is %g percent.\n', ...
    sqrt((mean(D_fit_Corr_2n3)-D_true)^2 + std(D_fit_Corr_2n3)^2)/D_true*100);
fprintf('The total number of trajectories is %g.\n', ns);
fprintf('The number of least square fit performed is %d.\n', num_group);
end

toc
%% save the results
% 
g_G_Corr_Corr3_Corr2n3_N2Tp10 = [D_fit_g2, D_fit_G, D_fit_Corr, D_fit_Corr3,...
    D_fit_Corr_2n3];
%save('D_fit.mat', 'g_G_Corr_Corr3_Corr2n3_N2T1')
%save('D_fit.mat', 'g_G_Corr_Corr3_Corr2n3_N2Tp01', '-append')
