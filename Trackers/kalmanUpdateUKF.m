function track = kalmanUpdateUKF(track, y, t, sensor_params)
    R = sensor_params.R;
    lim = sensor_params.lim;
    h = sensor_params.h;
    N = track.dimension;
    mod_index = sensor_params.modding_param_index;

    kappa = 3 - N;
    w = zeros(2*N + 1, 1);
    w(1) = kappa / (kappa + N);
    w(2:end) = 0.5 / (kappa + N);
    L = chol(track.P)';
    sigma_offset = sqrt(N + kappa) * L;
    pt_range = 1:N;
    x_sigmas = zeros(N, 2*N + 1);
    x_sigmas(:, 1) = track.x;
    x_sigmas(:, pt_range .* 2) = x_sigmas(:, 1) + sigma_offset;
    x_sigmas(:, pt_range .* 2 + 1) = x_sigmas(:, 1) - sigma_offset;

    y_sigmas = h(x_sigmas);
    y_bar = y_sigmas * w;
    % x_bar = x_sigmas * w;

    x_err = (x_sigmas - track.x);
    y_err = (y_sigmas - y_bar);
    Pxz = x_err * diag(w) * y_err';         % P*H'
    Pzz = y_err * diag(w) * y_err' + R;     % S
    K = Pxz / Pzz;

    innov = twoPiMod(y, y_bar, mod_index);

    track.x = track.x + K * innov;
    track.P = track.P - K * Pzz * K';
    track.P = (track.P + track.P')/2;

    track.t = t; 
    track = upgradeStatus(track, lim);

end