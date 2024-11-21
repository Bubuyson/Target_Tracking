function [track_id, is_new_track] = mapDist2Tracks(dists, track_ids, threshold)
    if isempty(dists) || min(dists) > threshold
        is_new_track = true;
        track_id = -1;
    else
        is_new_track = false;
        [~, id] = min(dists);
        track_id = track_ids(id);
    end  
    
end