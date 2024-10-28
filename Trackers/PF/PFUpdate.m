function track = PFUpdate(track, H, y, t, sensor_params)
    R = sensor_params.R;
    lim = sensor_params.lim;
    h = sensor_params.h;
    resample_method = sensor_params.resample_method;
    need_resample_method = sensor_params.need_resample_method;
    resample_parameter = sensor_params.resample_parameter;
    n_particle = length(track.mu);

    y_ = h(track.particles);
    S = H * track.P * H' + R;
    likelihoods = mvnpdf(y', y_', S)';
    track.mu = likelihoods .* track.mu;
    track.mu = track.mu ./ sum(track.mu);
    chosen_indices = 1:n_particle;                                      %#ok
    if need_resample_method == "MaxWeight" 
        if max(track.mu) > resample_parameter
            if resample_method == "Multinomial"
                cum_sum = cumsum(track.mu);
                cum_sum(end) = 1;
                rand_val = rand(1, n_particle);
                chosen_indices = zeros(1, n_particle);
                for i = 1:n_particle
                    chosen_indices(i) = binarySearchLess(cum_sum, rand_val(i));
                end
                track.mu = ones(1, n_particle) ./ n_particle;
            else
                error("This resample method is not implemented")
            end
        end
    else
        error("This need resample method is not implemented")
    end
    track.particles = track.particles(:, chosen_indices);
    % TODO: Modding of the states, one idea is to separately handle 
    % positive and negative values of azimuth, then combine
    track = combineParticles(track);
    track.t = t; 
    track = upgradeStatus(track, lim);


    function index = binarySearchLess(sortedArray, target)
        % Returns the index of the first element that is less than the searched
        % value, assuming target can not be more than the last element of the
        % array
        l = int32(1); r = int32(length(sortedArray));
        assert(target <= sortedArray(end), 'Choose target less than the maximum (end) value of the array')
        while l<=r
            m = (l + r) / 2;
            if(sortedArray(m) == target)
                index = m;
                break
            end
            if(sortedArray(m) > target)
                r = m-1;
                index = m;
            else
                l = m+1;
            end
        end
    end
end