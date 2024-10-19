function fa_meas = insertFA(meas, insertion, insert_indexes)
    [r, c] = size(meas);
    fa_meas = zeros(r, c + length(insert_indexes));
    itr = 1;
    prev_index = 1;
    if ~isempty(insert_indexes)
        next_index = insert_indexes(itr);
    else
        fa_meas = meas;
        return
    end
    while 1
        fa_meas(:, prev_index + itr - 1:next_index + itr - 1) = meas(:, prev_index:next_index);
        fa_meas(:, next_index + itr) = insertion(:, itr);
        itr = itr + 1;
        prev_index = next_index + 1;
        if itr <= size(insertion, 2)
            next_index = insert_indexes(itr);
        else
            next_index = size(meas, 2);
            fa_meas(:, prev_index + itr - 1:next_index + itr - 1) = meas(:, prev_index:next_index);
            break;
        end
           
    end

end