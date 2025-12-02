%% ========================================================================
%  Optimised Canny Edge Detection (OCED)Using a Hybrid of Ants and Bees
%  Colony Technique enhanced with Adaptive Histogram Equalisation for 
%  Feature Mapping and Extraction in Medical Images 
%  Author: Emmanuel, E.A.
%  ========================================================================
clear; clc; close all;
rng(0);
maxNumCompThreads(1);

%% ========================================================================
%  Folder Setup
%  ========================================================================
inputFolder = "inputFolder path";         % Replace with actual path
outputFolder = "OutputFolder name";       % Output folder for enhanced images

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% ========================================================================
%  Load Image Files
%  ========================================================================
imageFiles = [
    dir(fullfile(inputFolder, '*.jpg'))
    dir(fullfile(inputFolder, '*.jpeg'))
    dir(fullfile(inputFolder, '*.png'))
]; % You can add '*.tif','*.tiff', '*.bmp', etc.

numImages = min(length(imageFiles), 55);   % Limit to first 55 images

%% ========================================================================
%  Parameter Ranges
%  ========================================================================
cannyParamsRange = [0.010, 0.100; 0.101, 0.300]; 
sigmaRange      = [0.5, 1.5];
clipLimitRange  = [0.01, 1.00];
numTilesRange   = [4, 16];

%% ========================================================================
%  Optimization Setup
%  ========================================================================
numAnts = 10;
numBees = 10;
maxIterations = 50;

%% ========================================================================
%  Metrics Storage
%  ========================================================================
metricsAll = cell(numImages, 13);

%% ========================================================================
%  Main Processing Loop
%  ========================================================================
for idx = 1:numImages

    imagePath = fullfile(inputFolder, imageFiles(idx).name);
    originalImage = imread(imagePath);

    % Convert to grayscale
    if size(originalImage,3) == 3
        grayImage = rgb2gray(originalImage);
    else
        grayImage = originalImage;
    end

    % Traditional CED reference edge
    referenceEdge = edge(grayImage, 'canny');

    % Objective function for optimization
    objectiveFunction = @(params) evaluateParams(grayImage, params, referenceEdge);

    % Initialize ACO and BCO populations
    antSolutions = initializeSolutions(numAnts, cannyParamsRange, sigmaRange, clipLimitRange, numTilesRange);
    beeSolutions = initializeSolutions(numBees, cannyParamsRange, sigmaRange, clipLimitRange, numTilesRange);

    bestFitness = -inf;

    %% ====================================================================
    %  HABCO Hybrid Search
    %  ====================================================================
    for iter = 1:maxIterations

        % ACO phase
        for i = 1:numAnts
            newParams = perturbSolution(antSolutions(i,:), cannyParamsRange, sigmaRange, clipLimitRange, numTilesRange);
            fitness = objectiveFunction(newParams);

            if fitness > bestFitness
                bestFitness = fitness;
                bestSolution = newParams;
            end
        end

        % BCO phase
        for j = 1:numBees
            newParams = perturbSolution(beeSolutions(j,:), cannyParamsRange, sigmaRange, clipLimitRange, numTilesRange);
            fitness = objectiveFunction(newParams);

            if fitness > bestFitness
                bestFitness = fitness;
                bestSolution = newParams;
            end
        end
    end

    %% ====================================================================
    %  Final Optimised Edge Detection
    %  ====================================================================
    f = @() optimizedEdgeDetection(grayImage, bestSolution);
    executionTime = timeit(f);

    [~, enhancedEdges] = optimizedEdgeDetection(grayImage, bestSolution);

    %% ====================================================================
    %  Evaluation Metrics
    %  ====================================================================
    [precision, recall, f1Score, accuracy] = calculateScores(referenceEdge, enhancedEdges);

    mseValue     = immse(double(enhancedEdges), double(referenceEdge));
    mbeValue     = mean(double(enhancedEdges(:)) - double(referenceEdge(:)));
    psnrValue    = psnr(double(enhancedEdges), double(referenceEdge));
    entropyValue = entropy(enhancedEdges);
    kappaValue   = calculateKappaIndex(referenceEdge, enhancedEdges);
    [ssimValue, ~] = ssim(double(enhancedEdges), double(referenceEdge));

    preservationIndex = sum(referenceEdge(:) & enhancedEdges(:)) / sum(referenceEdge(:));

    %% ====================================================================
    %  Store Metrics
    %  ====================================================================
    metricsAll(idx,:) = {
        imageFiles(idx).name,...
        round(accuracy,4),...
        round(mseValue,4),...
        round(mbeValue,4),...
        round(psnrValue,4),...
        round(entropyValue,4),...
        round(kappaValue,4),...
        round(ssimValue,4),...
        round(precision,4),...
        round(recall,4),...
        round(f1Score,4),...
        round(preservationIndex,4),...
        round(executionTime,4)
    };

    % Save edge image
    [~, name, ~] = fileparts(imageFiles(idx).name);
    imwrite(enhancedEdges, fullfile(outputFolder, name + "_enhanced.png"));

end

%% ========================================================================
%  Save Metrics to Excel
%  ========================================================================
metricsTable = cell2table(metricsAll, ...
    'VariableNames', {'ImageName','Accuracy','MSE','MBE','PSNR_dB','Entropy', ...
                      'Kappa','SSIM','Precision','Recall','F1Score', ...
                      'PreservationIndex','ExecutionTime'});

writetable(metricsTable, 'metrics_evaluation.xlsx');

disp("Batch processing completed. Metrics saved.");

%% ========================================================================
%  Helper Functions
%  ========================================================================

function [edges, enhancedEdges] = optimizedEdgeDetection(grayImage, params)
    lowerThreshold = params(1);
    upperThreshold = params(2);
    sigma         = params(3);
    clipLimit     = params(4);
    numTiles      = round(params(5));

    edges = edge(grayImage, 'canny', [lowerThreshold, upperThreshold], sigma);

    enhancedEdges = adapthisteq(uint8(edges * 255), ...
        'ClipLimit', clipLimit, 'NumTiles', [numTiles numTiles]);

    enhancedEdges = imbinarize(enhancedEdges);
end

function [precision, recall, f1Score, accuracy] = calculateScores(reference, enhanced)
    TP = sum(reference(:) & enhanced(:));
    FP = sum(~reference(:) & enhanced(:));
    FN = sum(reference(:) & ~enhanced(:));
    TN = sum(~reference(:) & ~enhanced(:));

    precision = TP / (TP + FP + eps);
    recall    = TP / (TP + FN + eps);
    f1Score   = 2 * (precision * recall) / (precision + recall + eps);
    accuracy  = (TP + TN) / (TP + TN + FP + FN + eps);
end

function fitness = evaluateParams(image, params, referenceEdge)
    [~, enhancedEdges] = optimizedEdgeDetection(image, params);
    [~, ~, f1Score] = calculateScores(referenceEdge, enhancedEdges);
    fitness = f1Score;
end

function solutions = initializeSolutions(num, cRange, sRange, clRange, tRange)
    solutions = zeros(num,5);
    for i = 1:num
        solutions(i,:) = [
            randInRange(cRange(1,:)), ...
            randInRange(cRange(2,:)), ...
            randInRange(sRange), ...
            randInRange(clRange), ...
            randInRange(tRange)
        ];
    end
end

function newSolution = perturbSolution(solution, cRange, sRange, clRange, tRange)
    perturb = 0.1 * (rand(1,5) - 0.5);
    newSolution = solution + perturb;

    newSolution(1) = min(max(newSolution(1), cRange(1,1)), cRange(1,2));
    newSolution(2) = min(max(newSolution(2), cRange(2,1)), cRange(2,2));
    newSolution(3) = min(max(newSolution(3), sRange(1)), sRange(2));
    newSolution(4) = min(max(newSolution(4), clRange(1)), clRange(2));
    newSolution(5) = round(min(max(newSolution(5), tRange(1)), tRange(2)));
end

function val = randInRange(range)
    val = range(1) + (range(2) - range(1)) * rand;
end

function kappa = calculateKappaIndex(image1, image2)
    image1 = double(image1(:));
    image2 = double(image2(:));

    confusionMatrix = confusionmat(image1, image2);

    total = sum(confusionMatrix(:));
    po = sum(diag(confusionMatrix)) / total;
    pe = sum(sum(confusionMatrix,1).*sum(confusionMatrix,2)) / total^2;

    kappa = (po - pe) / (1 - pe + eps);
end
