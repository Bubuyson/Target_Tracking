function cart = sph2Cart_(sph, dims) 
    if nargin == 1
        dims = 1:3;
    end
    assert(size(sph, 1) == 3 || size(sph, 1) == 6, "Invalid dimesions")
    r = sph(1, :);
    phi = sph(2, :);
    theta = sph(3, :);
    x = r .* cos(theta) .* cos(phi);
    y = r .* cos(theta) .* sin(phi);
    z = r .* sin(theta);
    cart = [x; y; z];
    cart = cart(dims, :);
end

