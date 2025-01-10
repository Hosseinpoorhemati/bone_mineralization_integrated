function [peaks_locs, transition_point] = lag_time_calc(x, y)

    % Compute the derivative of y with respect to x
    dy_dx = diff(y) ./ diff(x);
    
    % Find peaks in the derivative to identify significant changes
    [peaks, locs] = findpeaks(dy_dx);
    if ~isempty(peaks)
        transition_point = x(locs(1) + 1); % Account for the diff operation
    else
        transition_point = NaN;
    end
    peaks_locs = [peaks, locs]

end
