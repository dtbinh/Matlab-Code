function varargout = histline( data , bins , plottype, lineFormat, smooth_WDW)
%This function plots a histogram of the data, but with a single line representing the data instead of bars 
%   It accepts the data, and the specific bins as input vectors
    
    if nargin < 2
        error('Please use at least two arguments');
    elseif ~(isvector(data) & isvector(bins))
        error('Data and Bins must both be vectors');
    end
    if nargin < 3
        plottype = 'none';
    end
    if nargin < 4
        lineFormat = 'b';
    end
    if nargin < 5
        smooth_WDW = 'dontSmooth';
    end
    data = squeeze(data); % Makes it easier
    bins=sort(bins); % Sorting bins into ascending order
    data(data<bins(1)) = NaN; % Removing anything outside the bins
    data_line=nan(length(bins),1);
    for n=2:length(bins)
        data_line(n-1) = length(find(data<bins(n)));
        data(data<bins(n)) = NaN;
    end
    if ~ischar(smooth_WDW) & smooth_WDW>0
        data_line=smooth(data_line,round(smooth_WDW));
    end
    if strcmp(plottype, 'stairs')
        handle = stairs(bins,data_line,lineFormat);
        varargout{1,2} = [data_line handle];
    elseif strcmp(plottype, 'plot')
        handle = plot(bins,data_line,lineFormat);
        varargout{1,2} = [data_line handle];
    elseif strcmp(plottype, 'none')
        varargout{1} = data_line;
    end
end

