%% pdfcounts
% Return data points for a PDF
%
% Tristan Abbott // Massachusetts Institute of Technology // 03/02/2017
%
%%% Syntax
%   [f, centers, edges, bin] = pdfcounts(X, ...)
%
%%% Description
% This function calculated a PDF based on some data by first constructing a
% histogram from the data (using the built-in histcounts function) and then
% normalizing the value of each histogram bin such that the values are
% probability densities. In other words, the PDF values and edges returned
% by the function are such that
%
%   sum((edges(2:end) - edges(1:end-1)).*f) == 1
%
%%% Inputs
% *X - data:* data from which to construct the PDF.
%
% *nbins (optional) - number of PDF points.*
%
% *edges (optional) - edges for PDF bins.*
%
%%% Outputs
% *f - PDF points:* PDF values calculated from the samples.
%
% *centers - PDF point locations:* approximate location for each PDF value.
% Calculated by averaging the left and right bin edge values.
%
% *edges - bin edges:* bin edge locations, either as specified, if the
% edges input argument is used, or as returned by histcounts.
%
% *bin - bin indices:* an array the same shape as X with values
% corresponding to the indices of the bins that elements of X were placed
% in.
%
% <../test/html/pdfcounts_test.html Tests>

function [f, centers, edges, bin] = pdfcounts(X, n)

    if nargin == 1
        [N, edges, bin] = histcounts(X);
    else
        [N, edges, bin] = histcounts(X, n);
    end
    w = edges(2:end) - edges(1:end-1);
    f = (N/numel(X))./w;
    centers = 0.5*(edges(2:end) + edges(1:end-1));

end