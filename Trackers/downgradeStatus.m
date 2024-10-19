function updated_track = downgradeStatus(track)
    updated_track = track;
    if track.status == "unconfirmed"
        updated_track.status = "to be deleted";
        assert(updated_track.confirmation_counter == 0, "There is a bug in downgradeStatus")
    elseif track.status== "tentative"
        updated_track.status = "unconfirmed";
        assert(updated_track.confirmation_counter == 0, "There is a bug in downgradeStatus")
    elseif track.status == "confirmed"
        if updated_track.confirmation_counter > 1
            updated_track.status = "confirmed";
            updated_track.confirmation_counter = updated_track.confirmation_counter - 1;
        else
            updated_track.status = "tentative";
            updated_track.confirmation_counter = updated_track.confirmation_counter - 1;
            assert(updated_track.confirmation_counter == 0, "There is a bug in downgradeStatus")
        end
    else
        error("Name error in downgradeStatus")
    end
end