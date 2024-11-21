function updated_track = upgradeStatus(track, confirmation_limit)
    updated_track = track;
    if track.status == "unconfirmed"
        updated_track.status = "tentative";
    elseif track.status == "tentative"
        updated_track.status = "confirmed";
        assert(updated_track.confirmation_counter == 0, "There is a bug in upgradeStatus")
        updated_track.confirmation_counter = 1;
    elseif track.status == "confirmed"
        updated_track.status = "confirmed";
        updated_track.confirmation_counter = min([updated_track.confirmation_counter + 1, confirmation_limit]);
    else
        error("Name error in upgradeStatus")
    end
end