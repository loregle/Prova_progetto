clearvars;close all;clc;
edges = [2.0 3.0 4.0 5.0 6.0]; % Define bin edges
counts = [2 4 3 1]; % Corresponding bin counts

% Reconstruct histogram
figure;
bar(edges(1:end-1), counts, 'histc'); % Use 'histc' to match bin edges
xlabel('Peak Values');
ylabel('Frequency');
title('Reconstructed Histogram');
