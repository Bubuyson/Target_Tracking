function out = poissrnd_NT(mu, num)
    arguments
        mu(1, 1) int64 = 5;
        num(1, :) int64 = 1;
    end
    out = -1*ones(1, num);
    for i = 1:num
        T = 0;
        while T <1
           E = -(1/double(mu))*log(rand);
           T = T+E;
           out(i) = out(i)+1;
        end
    end
end