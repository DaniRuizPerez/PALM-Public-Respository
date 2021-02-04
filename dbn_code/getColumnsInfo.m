function [cols_names, cols_types] = getColumnsInfo(cols, delimiter)
%
% Author: Jose Lugo-Martinez
% Advisor Profs: Giri Narasimhan and Ziv Bar-Joseph
% Description: Given a cell array of strings, where each string contains the 
% column attribute information: name and type. 
% Outpus each column attribute into a cell array of strings
%
% INPUT: 
%     colsAttributes % cell array of strings, where each string contains the 
%                    % column attribute information: name and type. 
%     delimiter      % user-specified delimiter. Default '$'. (optional)
%
% OUTPUT:
%     cols_names % vector of strings for column names
%     cols_types % vector of strings for column types
%
% (c) Jose Lugo-Martinez 2019.  MIT license. See cgbayesnets_license.txt.

    if (nargin < 2)
        delimiter = '$';
    end

    cols_names = string(zeros(1, length(cols))); %cell array version repmat({'0'}, 1, length(cols));
    cols_types = string(zeros(1, length(cols))); %cell array version repmat({'0'}, 1, length(cols));
    for i = (1:length(cols))   
        columnAttributes = strsplit(cols{i}, delimiter);
        cols_names(i) = columnAttributes{1};
        if length(columnAttributes) > 1
            cols_types(i) = columnAttributes{2};
        end
    end
end

