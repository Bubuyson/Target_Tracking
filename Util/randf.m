function val = randf(Arr, Len)
    arguments
        Arr(:, 2) double
        Len(1, 1) int64 = 1; 
    end
    val = (Arr(:,2) - Arr(:,1)).*rand(size(Arr, 1), Len) + Arr(:,1);
end