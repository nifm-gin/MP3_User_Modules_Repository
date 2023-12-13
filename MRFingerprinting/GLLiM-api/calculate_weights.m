function weights = calculate_weights(signal_value, bins, values)
    % Find the bin index where signal_value falls into
    bin_index = find(signal_value >= bins, 1, 'last');
    
    % If signal_value is greater than the largest bin, use the last value
    if isempty(bin_index)
        bin_index = length(bins);
    end
    
    % Get the weight based on the bin index
    weights = values(bin_index);
end