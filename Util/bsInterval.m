function [res, indexes]= bsInterval(arr, low, high)
    if low == high
        low_idx = binarySearchLower2(arr, low);
        high_idx = binarySearch(arr, high, "upper");
        if low_idx <= high_idx
            indexes = low_idx:high_idx;
            res = arr(indexes);
        else
            indexes = [];
            res = [];
        end
        return
    end
    low_idx = binarySearch(arr, low, "lower");
    high_idx = binarySearch(arr, high, "upper");
    indexes = low_idx:high_idx;
    res = arr(indexes);
    function idx = binarySearch(arr, value, bound)
        l = 1;
        r = length(arr);
        while l <= r
            m = floor((l + r) / 2);
            if strcmp(bound, "lower")
                if arr(m) > value
                    r = m - 1;
                else
                    l = m + 1;
                end
            
            elseif strcmp(bound, "upper")
                if arr(m) <= value
                    l = m + 1;
                else
                    r = m - 1;
                end
            end
        end
        if strcmp(bound, "upper")
            idx = r;
        elseif strcmp(bound, "lower")
            idx = l;
        end
    end

    function idx = binarySearchLower2(arr, value)
        l = 1;
        r = length(arr);
        while l <= r
            m = floor((l + r) / 2);
            if arr(m) >= value
                r = m - 1;
            else
                l = m + 1;
            end
        end
        idx = l;
    end

end