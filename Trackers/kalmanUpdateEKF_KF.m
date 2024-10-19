function track = kalmanUpdateEKF_KF(track, H, y, t, sensor_params)
    R = sensor_params.R;
    lim = sensor_params.lim;
    h = sensor_params.h;
    mod_index = sensor_params.modding_param_index;

    S = H * track.P * H' + R;
    K = track.P * H' * inv(S);    %#ok

    if (class(h) == "function_handle")
        innov = twoPiMod(y, h(track.x), mod_index);
    else
        innov = twoPiMod(y, H*track.x, mod_index);
    end
    track.x = track.x + K * innov;
    track.P = track.P - K * S * K';
    track.P = (track.P + track.P')/2;

    track.t = t; 
    track = upgradeStatus(track, lim);
end