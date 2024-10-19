function indexes = bsLess(arr, num)
    indexes = zeros(1, length(num));
    for i = 1:length(num)
        l = int32(1); r = int32(length(arr));
        indexes(i) = 0;
        while l<=r
            m = floor((l+r)/2);
            if arr(m) <= num(i)
                indexes(i) = m;
                l = m + 1;
            else
                r = m -1;
            end
        end    
    end

end