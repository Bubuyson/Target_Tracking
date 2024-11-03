function tracks = PFPredict(tracks, sensor_params, dt)
    assert(dt >= 0, "There is error with measurement sequence")
    if isempty(tracks) || dt == 0
        return
    end
    A = sensor_params.A(dt);
    B = sensor_params.B(dt);
    Q = sensor_params.Q;
    n_particles = length(tracks(1).mu);
    for i = 1:length(tracks)
        tracks(i).particles = A * tracks(i).particles + B * mvnrnd([0 0 0]', Q, n_particles)';
        tracks(i) = combineParticles(tracks(i));
    end
end