%% ========================================================================
%   MASTER SCRIPT FOR OCED SENSITIVITY AND ROBUSTNESS ANALYSIS
%   Author: Emmanuel E.A.
%   This script automatically runs:
%       1. Sensitivity analysis to HABCO parameters
%       2. Robustness analysis under different noise types
%       3. Generates heatmaps, convergence plots, robustness curves
%   ========================================================================

clc; clear; close all;
rng(0);
maxNumCompThreads(1);

%% ========================================================================
%  Load image
%  ========================================================================
imageFile = "C:\Users\User\Desktop\Medical Image data set\SampleXray\preprocessedimages_SampleXray\preprocessed_x-ray_8.jpeg";
img = imread(imageFile);

if size(img,3) == 3
    grayImage = rgb2gray(img);
else
    grayImage = img;
end

reference = edge(grayImage,"canny"); % Use as pseudo-ground truth

%% ========================================================================
%  SENSITIVITY ANALYSIS SETUP
%  ========================================================================
population_sizes = [5 10 20 30 40 50];
iteration_sizes   = [10 20 40 60 80 100];

F1_pop = zeros(length(population_sizes),1);
F1_iter = zeros(length(iteration_sizes),1);

fprintf("Running OCED Sensitivity Analysis...\n");

%% ========================================================================
%  Sweep 1: Sensitivity to Population Size
%  ========================================================================
for p = 1:length(population_sizes)
    [~, ~, f1] = OCED_Run(grayImage, reference, population_sizes(p), 50);
    F1_pop(p) = f1;
    fprintf("Population %d completed. F1 = %.4f\n", population_sizes(p), f1);
end

%% ========================================================================
%  Sweep 2: Sensitivity to Number of Iterations
%  ========================================================================
for k = 1:length(iteration_sizes)
    [~, ~, f1] = OCED_Run(grayImage, reference, 20, iteration_sizes(k));
    F1_iter(k) = f1;
    fprintf("Iterations %d completed. F1 = %.4f\n", iteration_sizes(k), f1);
end

%% ========================================================================
%  Robustness Analysis under Noise
%  ========================================================================

noise_levels = 0:0.02:0.20;
noise_types = {'gaussian','salt & pepper','speckle'};

Robustness = zeros(length(noise_levels), length(noise_types));

fprintf("\nRunning Noise Robustness Analysis...\n");

for n = 1:length(noise_levels)

    sigma = noise_levels(n);

    noisy_g = imnoise(grayImage, 'gaussian', 0, sigma);
    noisy_s = imnoise(grayImage, 'salt & pepper', sigma);
    noisy_k = imnoise(grayImage, 'speckle', sigma);

    imgs = {noisy_g, noisy_s, noisy_k};

    for t = 1:length(noise_types)
        [~, ~, f1] = OCED_Run(imgs{t}, reference, 20, 50);
        Robustness(n,t) = f1;
    end

    fprintf("Noise level %.2f completed\n", sigma);
end

%% ========================================================================
%  PLOTS AND VISUALISATIONS
%  ========================================================================

figure;
plot(population_sizes, F1_pop, '-o','LineWidth',1.5);
xlabel('Population Size'); ylabel('F1 Score');
title('Sensitivity of OCED to Population Size'); grid on;

figure;
plot(iteration_sizes, F1_iter, '-s','LineWidth',1.5);
xlabel('Number of Iterations'); ylabel('F1 Score');
title('Sensitivity of OCED to Iterations'); grid on;

% Robustness curves
figure; hold on;
plot(noise_levels, Robustness(:,1), '-o', 'LineWidth',1.5);
plot(noise_levels, Robustness(:,2), '-s', 'LineWidth',1.5);
plot(noise_levels, Robustness(:,3), '-^', 'LineWidth',1.5);
xlabel('Noise Level'); ylabel('F1 Score');
title('OCED Robustness under Noise');
legend('Gaussian','Salt & Pepper','Speckle');
grid on;

%% Heatmap for robustness
figure;
heatmap(noise_types, string(noise_levels), Robustness);
xlabel('Noise Type');
ylabel('Noise Level');
title('Heatmap of OCED Robustness');

fprintf("\nAll experiments completed successfully.\n");

