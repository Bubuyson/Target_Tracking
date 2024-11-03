function track = IMMUpdate(track, H, y, t, sensor_params)
    R = sensor_params.R;
    lim = sensor_params.lim;
    h = sensor_params.h;
    mod_index = sensor_params.modding_param_index;
    n_mode = track.n_mode;
    likelihood = zeros(1, n_mode);
    for i = 1:n_mode
        S = H * track.P(:, :, i) * H' + R;
        K = track.P(:, :, i) * H' * inv(S);    %#ok
        if (class(h) == "function_handle")
            y_ = h(track.x(:, i));
        else
            y_ = H*track.x(:, i);
        end
        innov = twoPiMod(y, y_, mod_index);
        track.x(:, i) = track.x(:, i) + K * innov;
        track.P(:, :, i) = track.P(:, :, i) - K * S * K';
        track.P(:, :, i) = (track.P(:, :, i) + track.P(:, :, i)')/2;

        likelihood(i) = mvnpdf(y, y_, S);
    end
    track.mu = likelihood .* track.mu / (sum(likelihood .* track.mu));
    track.t = t; 
    track = upgradeStatus(track, lim);
end