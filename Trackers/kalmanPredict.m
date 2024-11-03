function tracks = kalmanPredict(tracks, sensor_params, dt)
    assert(dt >= 0, "There is error with measurement sequence")
    if isempty(tracks) || dt == 0
        return
    end
    A = sensor_params.A(dt);
    B = sensor_params.B(dt);
    Q = sensor_params.Q;
    for i = 1:length(tracks)
        tracks(i).x = A * tracks(i).x;
        tracks(i).P = A * tracks(i).P * A' + B * Q * B';
    end
end