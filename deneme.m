clc; clear all; close all
% tracks = [1,1; 2,2];
% dets = [1.1, 1.1; 2.1, 2.1; 1.5, 3];
% for i = size(tracks, 1):-1:1
%     delta = dets - tracks(i, :);
%     costMatrix(i, :) = sqrt(sum(delta .^ 2, 2));
% end

costMatrix = [
            0.157969155795771   11.1071654829740
            12579099722307268  16          
            12  18.439438082324447
            14  14.4444444444444444444444444
            12.579099722307268  12.587724564818624];

% costMatrix = [
%             0.157969155795771   11.1071654829740
%             12579099722307268  16          ];

% costMatrix = [1.157969155795771   114.1071654829740
%                 12.579099722307268  12.587724564818624
%                 12  18.439438082324447
%                 11.9999999999999999  14.4444444444444444444444444
%                 12.579099722307268  12.587724564818624
%                 12.579099722307268  12.587724564818624
%                 12.579099722307268  12.587724564818624
%                 12.579099722307268  12.587724564818624
%                 12.579099722307268  12.587724564818624];
costofnonassignment = 1e10;
[assignments, unassignedTracks, unassignedDetections] = ...
    assignmunkres(costMatrix,costofnonassignment);

disp(assignments)
disp(unassignedTracks)
disp(unassignedDetections)



% t_arr = cell(1, n_sensor);
% data_arr = cell(1, n_sensor);
% for i = 1:t2t_hparams.n_fus
%     track = tracks{i};
%     for j = 1:length(track)
%         sensor_name = track{j}(1).sensor_name;
%         if sensor_name == vars.radar_name
%             t_arr{1}(end+1) = track{j}(end).t;
%             data_arr{1}(end+1) = length(track{j});
%         elseif sensor_name == vars.optic_name
%             t_arr{2}(end+1) = track{j}(end).t;
%             data_arr{2}(end+1) = length(track{j});
%         end
%     end
% end
% 
% figure(); hold on; grid on;
% for i = 1:n_sensor
%     sensor_name = track{i}(1).sensor_name;
%     if sensor_name == vars.radar_name
%         subplot(2, 1, 1); 
%         plot(t_arr{1}, data_arr{1}, 'LineWidth', 2); hold on; grid on
%         scatter(t_arr{1}, data_arr{1}, 15, 'filled', 'red');
%         title(vars.radar_name)
%     elseif sensor_name == vars.optic_name
%         subplot(2, 1, 2);
%         plot(t_arr{2}, data_arr{2}, 'LineWidth', 2); hold on; grid on
%         scatter(t_arr{2}, data_arr{2}, 15, 'filled', 'red');
%         title(vars.optic_name)
%     end
% end

% colors = chooseFromColorMap([], n_target);
% figure(); hold on; grid on; view(3)
% plot3(os.state(1, :), os.state(2, :), os.state(3, :), 'Color', os_color, 'LineWidth',2);
% for i = 1:n_target
% plot3(tgt_state{i}(1, :)...
% , tgt_state{i}(2, :) ...
% , tgt_state{i}(3, :)...
% ,'Color', colors(i, :), 'LineWidth', 2) % relative
% meas_cart_w_fa = sph2Cart_(meas.radar.sph{i});
% scatter3(meas_cart_w_fa(1, :) + os.state_meas_time(1).state(1, :)...
% , meas_cart_w_fa(2, :) + os.state_meas_time(1).state(2, :) ...
% , meas_cart_w_fa(3, :) + os.state_meas_time(1).state(3, :)...
% , 'MarkerEdgeColor', colors(i, :), 'Marker', '*') % relative
% end
% cubizer()
% radar_logs.confirmed_track_log 

% colors = chooseFromColorMap([], length(track_cell));
% figure(); hold on; grid on; view(3)
% plot3(os_state(1, :), os_state(2, :), os_state(3, :), 'Color', [0 0 0], 'LineWidth',2);
% for i = 1:n_target
%     plot3(tgt_state_meas_time{i}(1, :)...
%         , tgt_state_meas_time{i}(2, :) ...
%         , tgt_state_meas_time{i}(3, :)...
%         ,'Color', colors(i, :), 'LineWidth', 2)   % relative
%     scatter3(meas_cart_w_fa{i}(1, :) + os_state_meas_time_w_fa{i}(1, :)...
%         , meas_cart_w_fa{i}(2, :) + os_state_meas_time_w_fa{i}(2, :) ...
%         , meas_cart_w_fa{i}(3, :) + os_state_meas_time_w_fa{i}(3, :)...
%         , 'MarkerEdgeColor', colors(i, :), 'Marker', '*')   % relative
% end
% cubizer()
% radar_logs.confirmed_track_log
% 
