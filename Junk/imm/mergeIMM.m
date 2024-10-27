function [xhat,Phat] = mergeIMM(mu, x, P)
    Phat = zeros(size(x,1));
    xhat = x * mu'; 
    for j=1:length(mu)
        innov = x(:,j) - xhat;
        Phat = Phat + mu(j) * (P(:,:,j) + innov * innov');
    end

    function unitTest()
        k = [0.2 0.5 0.3];
        x = [5 3 6;
            4 12 7];
        P(:,:, 1) = eye(2);
        P(:,:, 2) = 10 * eye(2);
        P(:,:, 3) = 20 * eye(2);
        [xbar, Pbar] = mergeIMM(k, x, P);
    end    
end