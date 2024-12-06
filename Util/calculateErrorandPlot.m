function calculateErrorandPlot(radar_logs, meas, valid_n_monte, n_target, n_monte, algo_name, my_map)
    max_num_radar = 0;
    for m = 1:valid_n_monte
        if ~isempty(radar_logs{m})
            for i = 1:length(radar_logs{m}.confirmed_track_log)
                max_num_radar = max(max_num_radar, length(radar_logs{m}.confirmed_track_log{i}));
            end
        end
    end
    
    n_fusion = length(radar_logs{1}.confirmed_track_log); % Assumes all measurements are taken at the same time
    track_times = zeros(1, n_fusion); 
    for i = 1:n_fusion
        if ~isempty(radar_logs{m}.confirmed_track_log{i})
            track_times(i) = max([radar_logs{m}.confirmed_track_log{i}.t]);
        end
    end
    
        
    cart_error_prism = 0 * ones(n_target, n_fusion, n_monte);
    vel_error_prism = 0 * ones(n_target, n_fusion, n_monte);
    n_track_log = zeros(1, n_fusion);
    for m = 1:valid_n_monte
        for i = 1:n_fusion  
            n_track = length(radar_logs{m}.confirmed_track_log{i});
            n_track_log(i) = n_track;
            if n_track == 0
                continue
            end
            dists = zeros(n_target, n_track);
            diff_vectors = cell(1, n_target);
            for j = 1:n_target
                real_tgt_state = meas(m).radar.true_cart_all_meas_time{j}(:, i);
                diff_vectors{j} = zeros(6, n_track);
                for k = 1:n_track
                    diff_vectors{j}(:, k) = (radar_logs{m}.confirmed_track_log{i}(:, k).x - real_tgt_state);
                end
                dists(j, :) = sqrt(sum(diff_vectors{j}(1:3, :).^2));
            end
            [sensor_to_system, ~, ~] = assign(dists, 1e20);


            for k = 1:min(n_target, n_track)
               pos_dif = diff_vectors{sensor_to_system(k, 1)}(1:3, sensor_to_system(k, 2));
               vel_dif = diff_vectors{sensor_to_system(k, 1)}(4:6, sensor_to_system(k, 2));

               % range_error = sqrt(sum(pos_dif.^2));
               cartesian_pos_error = sqrt(sum(pos_dif.^2));
               cart_error_prism(sensor_to_system(k, 1), i, m) = cartesian_pos_error;

               % range_error = sqrt(sum(pos_dif.^2));
               cartesian_vel_error = sqrt(sum(vel_dif.^2));
               vel_error_prism(sensor_to_system(k, 1), i, m) = cartesian_vel_error;
            end
        end
    end
    
    cart_pos_error_mean_of_means = mean(sum(cart_error_prism, 1) ./ n_track_log, 3);
    cart_vel_error_mean_of_means = mean(sum(vel_error_prism, 1) ./ n_track_log, 3);
    non_zero_indexes = ~isnan(cart_pos_error_mean_of_means);
    
    
    cart_pos_error_mean_of_means = cart_pos_error_mean_of_means(non_zero_indexes);
    cart_vel_error_mean_of_means = cart_vel_error_mean_of_means(non_zero_indexes);
    track_times = track_times(non_zero_indexes);
    
    mean_cart_pos_error = mean(cart_pos_error_mean_of_means);
    mean_cart_vel_error = mean(cart_vel_error_mean_of_means);

    disp(my_map(algo_name) + " algorithm makes " + num2str(mean_cart_pos_error) + " m mean positional error")
    disp(my_map(algo_name) + " algorithm makes " + num2str(mean_cart_vel_error) + " m/s mean velocity error")
    
    figure();
    subplot(2, 1, 1); hold on; grid on
    plot(track_times, cart_pos_error_mean_of_means, 'LineWidth', 2, 'Color', 'r')
    xlabel('time (s)'); ylabel('Position error (m)'); title('Cartesian Position Errors wrt Time')
    subplot(2, 1, 2); hold on; grid on
    plot(track_times , cart_vel_error_mean_of_means, 'LineWidth', 2, 'Color', 'r')
    xlabel('time (s)'); ylabel('Velocity error (m/s)'); title('Cartesian Velocity Errors wrt Time')
    str = (my_map(algo_name) + " Results: pos error = " + num2str(mean_cart_pos_error) + ", vel error = " +  num2str(mean_cart_vel_error));
    sg = sgtitle(str);
    sg.FontSize = 12;

end