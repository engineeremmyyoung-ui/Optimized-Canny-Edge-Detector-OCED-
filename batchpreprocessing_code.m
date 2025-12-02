% =========================================================================
% Batch Preprocessing of Medical Images
% Steps: Grayscale Conversion, Resize, Normalization, Median Denoising
%
% Author: Emmanuel, E.A.
% Date: 5th December, 2025
% Description:
%   This script loads medical images from an input directory and applies
%   standard preprocessing steps suitable for medical image analysis tasks.
%
%   Outputs are saved in the specified output directory using the prefix
%   "preprocessed_".
%
% Usage:
%   1. Set inputFolder and outputFolder.
%   2. Run the script.
%
% Supported Formats: JPG, JPEG, PNG, BMP, TIFF
% =========================================================================

clc; clear; close all;

%% User Settings

inputFolder  = "INPUT_FOLDER_PATH";     % <-- Replace with your folder
outputFolder = "OUTPUT_FOLDER_PATH";    % <-- Replace with your folder

%% Validate Directories

if ~isfolder(inputFolder)
    error("Input folder not found: %s", inputFolder);
end

if ~isfolder(outputFolder)
    mkdir(outputFolder);
    fprintf("Created output folder: %s\n", outputFolder);
end

%% Supported Image Extensions

extensions = {"*.jpg", "*.jpeg", "*.png", "*.bmp", "*.tif", "*.tiff"};
imageFiles = [];

for k = 1:length(extensions)
    imageFiles = [imageFiles; dir(fullfile(inputFolder, extensions{k}))]; %#ok<AGROW>
end

if isempty(imageFiles)
    error("No supported image files found in the input folder.");
end

fprintf("Found %d images to preprocess...\n", length(imageFiles));

%% Processing Loop

for i = 1:length(imageFiles)

    filename = imageFiles(i).name;
    filepath = fullfile(inputFolder, filename);

    % Read image
    img = imread(filepath);

    % Convert RGB to grayscale if needed
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    % Resize image to 256x256
    resizedImg = imresize(img, [256 256]);

    % Normalize to [0, 1]
    normImg = mat2gray(resizedImg);

    % Noise reduction (median filter)
    denoisedImg = medfilt2(normImg, [3 3]);

    % Save output
    outputName = fullfile(outputFolder, ['preprocessed_' filename]);
    imwrite(denoisedImg, outputName);

    % Console feedback
    fprintf("Processed (%d/%d): %s\n", i, length(imageFiles), filename);

end

fprintf("\nBatch preprocessing complete.\n");
