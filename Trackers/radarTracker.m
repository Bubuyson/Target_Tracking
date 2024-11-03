function [radar_logs] = radarTracker(radar_meas, radar_hparams, vars)
    prev_t = 0;
    track_counter = 0;
    n_meas = size(radar_meas.sph_all_w_fa, 2);
    radar_logs.track_log = cell(1, n_meas);
    radar_logs.dists_log = cell(1, n_meas);
    radar_logs.confirmed_track_log = cell(1, n_meas);
    for i = 1:n_meas
        if ~exist('radar_prev_tracks', 'var')
            radar_prev_tracks = sensorTrackConstuctor();
        end
        y = radar_meas.sph_all_w_fa(:, i);
        t = radar_meas.all_times_w_fa(i);
        dt = t - prev_t;
        meas_cart = sph2Cart_(y, 1:3);
    
        radar_tracks = deleteTracks(radar_prev_tracks, t, radar_hparams.track_delete_time); %% includes downgrading
        clear radar_prev_tracks 
        radar_tracks = kalmanPredict(radar_tracks, radar_hparams, dt);  
        H = getNumericalJacobian(radar_hparams.h, meas_cart, [1e-4, 1e-4, 1e-4]');
        H_extended = [H, zeros(3)];
        [dists, dist_ids] = checkDist2Tracks(radar_tracks, y, radar_hparams, H_extended);
        [track_id, is_new_track] = mapDist2Tracks(dists, dist_ids, radar_hparams.track_threshold);
        if is_new_track
            inv_H = getNumericalJacobian(radar_hparams.inv_h, y, [1e-4, 1e-5, 1e-5]');
            track_counter = track_counter + 1;
            P = blkdiag(inv_H * radar_hparams.R * inv_H', (radar_hparams.v_max / radar_hparams.kappa)^2 * eye(3));
            track = sensorTrackConstuctor('sensor_name', vars.radar_name,...
                        'track_id', track_counter, ...
                         'sensor_id', vars.radar_id,...
                         'x', [meas_cart; radar_hparams.init_vel], 'P', P, ...
                         't', t, 'status', 'tentative', 'confirmation_counter',...
                         0, 'state_names', ["x pos", "y pos", "z pos",...
                         "x vel", "y vel", "z vel",], 'dimension', 6,...
                         "algo_name", radar_hparams.algo);
            radar_tracks(end+1) = track;      %#ok
        else
            for j = 1:length(radar_tracks)
                if radar_tracks(j).track_id == track_id
                    if radar_hparams.algo == vars.ekf_name
                        radar_tracks(j) = kalmanUpdateEKF_KF(radar_tracks(j), H_extended, y, t, radar_hparams);
                    elseif radar_hparams.algo == vars.ukf_name
                        radar_tracks(j) = kalmanUpdateUKF(radar_tracks(j), y, t, radar_hparams);
                    else
                        error("Update algorithm is not named correctly")
                    end
                    break
                end
            end
        end
        radar_logs.track_log{i} = radar_tracks;
    
        confirmed_tracks = sensorTrackConstuctor();
        for j = 1:length(radar_tracks)
            if radar_tracks(j).status == "confirmed" && radar_tracks(j).confirmation_counter >= radar_hparams.confirmation_number_req
                confirmed_tracks(end+1) = radar_tracks(j);            %#ok
            end
        end
        radar_logs.confirmed_track_log{i} = confirmed_tracks;
    
        radar_logs.dists_log{i} = dists;
        prev_t = t;
        radar_prev_tracks = radar_tracks;
        clear radar_tracks
    end

end