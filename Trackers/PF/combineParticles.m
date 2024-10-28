function track = combineParticles(track)
    n_particles = length(track.mu);
    track.x = track.particles * track.mu';
    track.P = zeros(size(track.P)); 
    for j = 1:n_particles
        innov = track.particles(j) - track.x;
        track.P = track.P + innov * innov';
    end
    track.P = track.P / n_particles;
end