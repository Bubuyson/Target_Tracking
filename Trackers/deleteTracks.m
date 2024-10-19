function updated_prev_tracks = deleteTracks(prev_tracks, t, track_delete_time)
    updated_prev_tracks = sensorTrackConstuctor();
    for i = 1:length(prev_tracks)
        if ~ ((t - prev_tracks(i).t) > track_delete_time)
            updated_prev_tracks(end+1) = prev_tracks(i);           %#ok
        else
            prev_tracks(i) = downgradeStatus(prev_tracks(i));
            if ~(strcmp(prev_tracks(i).status, "to be deleted"))
                updated_prev_tracks(end+1) = prev_tracks(i);           %#ok
            end
        end
    end
end