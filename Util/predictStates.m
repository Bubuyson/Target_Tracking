function tracks = predictStates(tracks, sensor_params, dt)
    A = sensor_params.A(dt);
    B = sensor_params.B(dt);
    Q = sensor_params.Q;
    for i = 1:length(tracks)
        tracks(i).x = A * tracks(i).x;
        tracks(i).P = A * tracks(i).P * A' + B * Q * B';
    end
end