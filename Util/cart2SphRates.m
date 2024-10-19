function sph = cart2SphRates(cart)
    assert(size(cart, 1) == 6, "Invalid dimesions")
    x  = cart(1, :);
    y  = cart(2, :);
    z  = cart(3, :);
    vx = cart(4, :);
    vy = cart(5, :);
    vz = cart(6, :);

    r_gnd = sqrt(x.^2 + y.^2);
    r     = sqrt(x.^2 + y.^2 + z.^2);

    phi = atan2(y, x);
    theta = atan2(z, r_gnd);

    vphi = zeros(size(x));
    vtheta = zeros(size(x));
    non_zero_gnd = (r_gnd ~= 0);
    vphi(non_zero_gnd) =  (x(non_zero_gnd) .* vy(non_zero_gnd) ...
        -y(non_zero_gnd) .* vx(non_zero_gnd)) ...
        ./ (r_gnd(non_zero_gnd).^2);

    non_zero_r = (r ~= 0);

    % vtheta(non_zero_r) = (vz(non_zero_r) ./ r(non_zero_r)) - ...
    %     (z(non_zero_r) .* (x(non_zero_r) .* vx(non_zero_r) + y(non_zero_r) .* vy(non_zero_r) ...
    %     + z(non_zero_r) .* vz(non_zero_r))) ./ r(non_zero_r).^2;
    vtheta(non_zero_r) = (x(non_zero_r) .* vy((non_zero_r)) - y(non_zero_r) .* vx(non_zero_r)) ./ (r_gnd(non_zero_r).^2 + z.^2);

    at_origin = (r==0);
    phi(at_origin) = 0;
    theta(at_origin) = 0;
    vphi(at_origin) = 0;
    vtheta(at_origin) = 0;

    sph = [phi; theta; vphi; vtheta];

    function sph = makeUnitTests()          %#ok
        carts = zeros(6, 1);                % -> [0 0 0 0]'                 origin, stationary
        carts = [1 0 0 1 0 0]';             % -> [0 0 0 0]'                 x axis pt go along x axis
        carts = [0 1 0 0 1 0]';             % -> [pi/2 0 0 0]'              y axis pt go along y axis
        carts = [0 0 1 0 0 1]';             % -> [0 pi/2 0 0]'              z axis pt go along z axis
        carts = [1 1 1 1 1 1]';             % -> [pi/4 0.6155 0 0]'         general point
        carts = [1 0 0 0 0 1]';             % -> [0 0 0 1]'                 x axis pt and going along z axis

        cart2SphRates(carts)
    end
end


    % vtheta(non_zero_r) = (r_gnd(non_zero_r) .* vz(non_zero_r) - z(non_zero_r) ...
    %     .* (x(non_zero_r) .* vx(non_zero_r) + y(non_zero_r) .* vy(non_zero_r)))...
    %     ./ r(non_zero_r).^2;