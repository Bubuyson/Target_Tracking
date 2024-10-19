function arr = cell2Mat_(cell_arr)
    arr = [];
    for i = 1:length(cell_arr)
        arr = [arr cell_arr{i}];         %#ok
    end

end