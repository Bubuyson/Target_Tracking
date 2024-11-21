function [radar_logs] = radarTrackerIMM(radar_meas, radar_hparams, vars)
    %% This function assumes that only Q values of the EKF trackers are different 
    prev_t = 0;
    track_counter = 0;
    n_meas = size(radar_meas.sph_all_w_fa, 2);
    radar_logs.track_log = cell(1, n_meas);
    radar_logs.dists_log = cell(1, n_meas);
    radar_logs.confirmed_track_log = cell(1, n_meas);
    for i = 1:n_meas
        if ~exist('radar_prev_tracks', 'var')
            radar_prev_tracks = sensorTrackConstuctorIMM();
        end
        y = radar_meas.sph_all_w_fa(:, i);
        t = radar_meas.all_times_w_fa(i);
        dt = t - prev_t;
        meas_cart = sph2Cart_(y, 1:3);
    
        radar_tracks = deleteTracksIMM(radar_prev_tracks, t, radar_hparams.track_delete_time); %% includes downgrading
        clear radar_prev_tracks 
        radar_tracks = IMMPredict(radar_tracks, radar_hparams, dt);  
        H = getNumericalJacobian(radar_hparams.h, meas_cart, [1e-4, 1e-4, 1e-4]');
        H_extended = [H, zeros(3)];
        %% TODO: Change this to mixing, this is not correct
        merged_tracks = mergeTracks(radar_tracks);
        [dists, dist_ids] = checkDist2Tracks(merged_tracks, y, radar_hparams, H_extended);
        [track_id, is_new_track] = mapDist2Tracks(dists, dist_ids, radar_hparams.track_threshold);
        if is_new_track
            inv_H = getNumericalJacobian(radar_hparams.inv_h, y, [1e-4, 1e-5, 1e-5]');
            track_counter = track_counter + 1;
            P = blkdiag(inv_H * radar_hparams.R * inv_H', (radar_hparams.v_max / radar_hparams.kappa)^2 * eye(3));
            track = sensorTrackConstuctorIMM('sensor_name', vars.radar_name,...
                        'track_id', track_counter, ...
                         'sensor_id', vars.radar_id,...
                         'n_mode', radar_hparams.n_mode, 'x', repmat([meas_cart; radar_hparams.init_vel], 1, 3), ...
                         'mu', ones(1, radar_hparams.n_mode) / radar_hparams.n_mode, 'P', repmat(P, 1, 1, 3), ...
                         't', t, 'status', 'tentative', 'confirmation_counter',...
                         0, 'state_names', ["x pos", "y pos", "z pos",...
                         "x vel", "y vel", "z vel",], 'dimension', 6,...
                         "algo_name", radar_hparams.algo);
            radar_tracks(end+1) = track;      %#ok
        else
            for j = 1:length(radar_tracks)
                if radar_tracks(j).track_id == track_id
                    radar_tracks(j) = IMMUpdate(radar_tracks(j), H_extended, y, t, radar_hparams);
                    break
                end
            end
        end
        radar_logs.track_log{i} = radar_tracks;
    
        confirmed_tracks = struct('track_id', {}, 't', {}, 'x', {}, 'mu', {}, 'P', {});
        for j = 1:length(radar_tracks)
            if radar_tracks(j).status == "confirmed" && radar_tracks(j).confirmation_counter >= radar_hparams.confirmation_number_req
                n = length(confirmed_tracks);
                merged_track = mergeTracks(radar_tracks(j));
                assert(~any(isnan(radar_tracks(j).x), 'all'))
                assert(~any(isnan(radar_tracks(j).mu), 'all'))
                confirmed_tracks(n + 1).x = merged_track.x;              
                confirmed_tracks(n + 1).P = merged_track.P;
                confirmed_tracks(n + 1).t = radar_tracks(j).t;            
                confirmed_tracks(n + 1).track_id = radar_tracks(j).track_id; 
                confirmed_tracks(n + 1).mu = radar_tracks(j).mu ; 
            end
        end
        radar_logs.confirmed_track_log{i} = confirmed_tracks;
    
        radar_logs.dists_log{i} = dists;
        prev_t = t;
        radar_prev_tracks = radar_tracks;
        clear radar_tracks
    end

end