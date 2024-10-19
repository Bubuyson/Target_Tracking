function out = mvnrnd_NT(mu, sigma, num)
    arguments
        mu(1, :) double = [0 0 0]
        sigma(:, :) double = [1 0 0; 0 1 0; 0 0 1]
        num(1, 1) int32 = 1
    end
    R = sqrt(sigma);
    out = repmat(mu,num,1) + randn(num, size(mu, 2))*R;
end
