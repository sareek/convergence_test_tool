close all; clc; clear all;
%chnage code in line 28 to 30 for selecting appropriate distribution, also
%change line 78 and 79
%change 
T1 = 100; T2 = 1000; T3 = 10000;

%Parameters
%a
mu = 2; sigma = sqrt(2);
%b
a = 2; b = 4;
%c
lambda = 2;

%%MATLAB inbuilt function for generating normal RV
norm_T1 = normrnd(2,sqrt(2), T1, 1);
norm_T2 = normrnd(2, sqrt(2), T2, 1);
norm_T3 = normrnd(2, sqrt(2), T3, 1);

% Rejection method
norm_f = @(x) (1/(sqrt(2*pi)*sigma))*exp(-((x-mu).^2)/(2*sigma^2));
%g function chosen as uniform function
norm_g = @(x) (x >= -5 & x <= 5) * (1/10); 
grnd_norm = @() -5 + 10 * rand(); 

%function reference: https://www.mathworks.com/help/stats/generating-random-data.html#br5k9hi-4
%comment/uncomment these for T1, T2, and T3 respectively
% X_norm = rejection_fun(norm_f, norm_g, grnd_norm, T1, 1);
% X_norm = rejection_fun(norm_f, norm_g, grnd_norm, T2, 1);
X_norm = rejection_fun(norm_f, norm_g, grnd_norm, T3, 1);

%histogram for normal variables with mean 2 and variance 2
% % figure(1);
% % h1 = histogram(norm_T1);
% % h1.BinWidth = 0.25;
% % hold on;

%%%%%Uniform random variable %%%%%%%%%%%

%MATLAB inbuilt funtion
unif_T1 = unifrnd(a, b, T1, 1);
unif_T2 = unifrnd(a, b, T2, 1);
unif_T3 = unifrnd(a, b, T3, 1);

% Rejection method

unif_f = @(x) (x >= a & x <= b) * (1 / (b - a));
%g function chosen as exponential function
unif_g = @(x) exp(-x);
grnd_unif = @() exprnd(1); 
% X_unif = rejection_fun(unif_f, unif_g, grnd_unif, T1, 1);
% X_unif = rejection_fun(unif_f, unif_g, grnd_unif, T2, 1);
X_unif = rejection_fun(unif_f, unif_g, grnd_unif, T3, 1);

%%%%%Exponential random variable %%%%%%%%%%%

% MATLAB inbuilt
exp_T1 = exprnd(1/2, T1, 1);
exp_T2 = exprnd(1/2, T2, 1);
exp_T3 = exprnd(1/2, T3, 1);

%Rejection
%exponential pdf
exp_f = @(x) lambda * exp(-lambda * x) .* (x >= 0); 
%g function chosen as normal function
exp_g= @(x) (1/(sqrt(2*pi))) * exp(-(x.^2) / 2); 
grnd_exp = @() randn();  
% X_exp = rejection_fun(exp_f, exp_g, grnd_exp, T1, 1);
% X_exp = rejection_fun(exp_f, exp_g, grnd_exp, T2, 1);
X_exp = rejection_fun(exp_f, exp_g, grnd_exp, T3, 1);


% Figures
type = {'Normal','Uniform','Exponential'};
bin_sizes = [3, 1, 0.5, 0.25];
data_rej = {X_norm, X_unif, X_exp};
% data_matlab = {norm_T1, unif_T1, exp_T1};
% data_matlab = {norm_T2, unif_T2, exp_T2};
data_matlab = {norm_T3, unif_T3, exp_T3};

for i = 1:length(type)
    figure;
    for j = 1:length(bin_sizes)
        subplot(2, length(bin_sizes), j);
        histogram(data_rej{i}, 'BinWidth', bin_sizes(j));
        title([type{i}, ' Bin Size=', num2str(bin_sizes(j))]);
        
        subplot(2, length(bin_sizes), j + length(bin_sizes));
        histogram(data_matlab{i}, 'BinWidth', bin_sizes(j));
        title([type{i}, ' Bin Size=', num2str(bin_sizes(j))]);
    end
end

disp('Estimation of the parameters of each of the population for each of the Tjs');
actual_params = {[mu, sigma^2], [a, b], [1/lambda]};
for i = 1:length(type)
    % estimate for all a,b,c
    if i == 1 
        estimated_mu = mean(data_rej{i});
        estimated_sigma = var(data_rej{i});
        estimated_params = [estimated_mu, estimated_sigma];
    elseif i == 2 
        estimated_a = min(data_rej{i});
        estimated_b = max(data_rej{i});
        estimated_params = [estimated_a, estimated_b];
    elseif i == 3 
        estimated_lambda = 1 / mean(data_rej{i});
        estimated_params = [1 / estimated_lambda];
    end
    
    
    % Compare to theoretical parameters
    disp(['Distribution: ', type{i}]);
    disp(['Actual Parameters (based on mathematical distribution): ', mat2str(actual_params{i})]);
    disp(['Estimated Parameters: ', mat2str(estimated_params)]);
    disp('  ');
end


function X = rejection_fun(f, g, generator_fun, row, col)
    X = zeros(row, col);
    for i = 1:row * col
        accept = false; 
        while ~accept
            u = rand(); % random variable
            v = generator_fun(); % rand variable from fun (g)
            %check if condition satisfies, if yes set accept true 
            if u <= f(v) / g(v)
                X(i) = v;
                accept = true;
            end
        end
    end
end
