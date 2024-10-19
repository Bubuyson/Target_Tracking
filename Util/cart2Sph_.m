function sph = cart2Sph_(cart, dims)
    if nargin == 1
        dims = 1:3;
    end
    assert(size(cart, 1) == 3 || size(cart, 1) == 6, "Invalid dimesions")
    x = cart(1, :);
    y = cart(2, :);
    z = cart(3, :);
    r = sqrt(x.^2 + y.^2 + z.^2);
    phi = atan2(y, x);
    theta = atan2(z, sqrt(x.^2 + y.^2));
    sph = [r; phi;theta];
    sph = sph(dims, :);
end


