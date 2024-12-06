function tracks = IMMPredict(tracks, sensor_params, dt)
    assert(dt >= 0, "There is error with measurement sequence")
    if isempty(tracks) || dt == 0
        return
    end
    A = sensor_params.A(dt);
    B = sensor_params.B(dt);
    for i = 1:length(tracks)
        for j = 1:sensor_params.n_mode
            Q = sensor_params.Q{j};
            tracks(i).x(:, j)  = A * tracks(i).x(:, j);
            tracks(i).P(:, :, j)  = A * tracks(i).P(:, :, j) * A' + B * Q * B';
        end
        tracks(i).mu = ((sensor_params.PI)^dt * tracks(i).mu')';
    end
end
