function [dists, ids] = checkDist2Tracks(tracks, y, sensor_params, H)
    dists = zeros(1, length(tracks));
    ids = zeros(1, length(tracks));
    h = sensor_params.h;
    mod_index = sensor_params.modding_param_index;
    R = sensor_params.R;
    for i = 1:length(tracks)
        if (class(h) == "function_handle")
            innov = twoPiMod(y, h(tracks(i).x), mod_index);
        else
            innov = twoPiMod(y, H*tracks(i).x, mod_index);
        end
        S = H * tracks(i).P * H' + R;
        dists(i) = sqrt(innov' / S * innov);  % no need for sqrt
        ids(i) = tracks(i).track_id;
        clear innov S
    end

end