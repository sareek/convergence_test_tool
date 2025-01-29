close all; clc; clear all;
Tj = [100, 1000, 10000]; num_bins = 30; 

% Generate distributions and compute sample means
for i = Tj
    % Random variables from part 1, normal, uniform and exponential
    %part 1
    norm_X = normrnd(2, sqrt(2), i, 1);
    unif_X = unifrnd(2, 4, i, 1);
    exp_X = exprnd(1/2, i, 1);

    % get normalized sum to get Yj for all three kinds of distribution
    norm_Yj = mean(norm_X);
    unif_Yj = mean(unif_X);
    exp_Yj = mean(exp_X);

    figure;
    subplot(3, 1, 1);
    %histogram for normalized sum of normal Xi, for i number of samples
    histogram(mean(normrnd(2, sqrt(2), [i, i])), num_bins);
    title(['Yj Normal for i = ', num2str(i)]);
    xlabel('Yj');
    ylabel('pdf');

    subplot(3, 1, 2);
    %histogram for normalized sum of uniform Xi, for i number of samples
    histogram(mean(unifrnd(2, 4, [i, i])), num_bins);
    title(['Yj Uniform for i = ', num2str(i)]);
    xlabel('Yj');
    ylabel('pdf');

    subplot(3, 1, 3);
    %histogram for normalized sum of exponential Xi, for i number of samples
    histogram(mean(exprnd(1/2, [i, i])), num_bins);
    title(['Yj Exponential for i = ', num2str(i)]);
    xlabel('Yj');
    ylabel('pdf');
end
