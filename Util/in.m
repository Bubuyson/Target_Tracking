function flag = in(arr, val)
    flag = ~isempty(find(arr == val, 1));
end
