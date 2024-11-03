function merged_tracks = mergeTracks(radar_tracks)
    if isempty(radar_tracks)
        merged_tracks = radar_tracks;
        return
    end
    n_track = length(radar_tracks);
    merged_tracks = mergedTrackConstuctorIMM();
    for i = 1:n_track
        merged_tracks(i).x = radar_tracks(i).x * radar_tracks(i).mu';
        merged_tracks(i).track_id = radar_tracks(i).track_id;
        merged_tracks(i).P = zeros(radar_tracks(i).dimension);
        for j = 1:radar_tracks(i).n_mode
            innov = merged_tracks(i).x - radar_tracks(i).x(:, j);
            merged_tracks(i).P = merged_tracks(i).P + radar_tracks(i).mu(j) ...
            * (radar_tracks(i).P(:, :, j) + innov * innov');
        end
    end
end