clc; clear; close all
%% Paths
addpath(genpath('Util'))
addpath(genpath('Trackers'))

%% Global variables
vars.radar_name  = "Radar" ;
vars.radar_id    = 1;
vars.ekf_name    = "EKF";
vars.ukf_name    = "UKF";
vars.imm_name    = "IMM";
vars.pf_name     = "PF";

%% Seed parameters
time_seed = 5;
is_seeded = true;
if is_seeded
    other_seed = 10;
end

%% Simulation parameters
n_monte = 10;
n_target = 5;
t_max = 50;
sim_dt = 0.01;
sim_t_arr = 0:sim_dt:t_max;
if is_seeded
    rng(other_seed)
end

%% Sensor Choices
sensor_options.radar = true;            % Field names are sensor specific and case sensitive

%% Mc run that will be plotted in the figure
index_mc = 1;
assert(index_mc <= n_monte, "Index of the mc that will be plotted should be less than or equal to the number of mc runs")

%% Motion models and noise matrix
const_vel3 = @(dt) [eye(3) eye(3) * dt; zeros(3) eye(3)];
global_motion_model = const_vel3(sim_dt);

B3 = @(dt) [eye(3) * dt^2/2; eye(3)* dt];

%%
sensor_list = createSensorList(sensor_options,vars);
include_radar = in(sensor_list, vars.radar_name);
n_sensor = length(sensor_list);

%% Radar hyperparameters
if include_radar
    %% Measurement time parameters
    radar_hparams.t_mean = 3;
    radar_hparams.t_var = 1;
    radar_hparams.bias = [0, 0, 0];         % [range az el]

    %% Track parameters
    % radar_hparams.track_threshold = 5e5;            % Beta
    radar_hparams.track_threshold = 25e0;            % Beta
    radar_hparams.track_delete_time = 5;            % Delete after x seconds of no updating
    radar_hparams.lim = 1;                         % confirmation limit
    radar_hparams.confirmation_number_req = 1;      % Required number of updates for confirmation

    %% False alarm parameters
    radar_hparams.lambda = 10;

    %% EKF UKF IMM selection        params = ["EKF", "UKF", "PF", "IMM"]
    radar_hparams.algo = vars.imm_name;

    %% Process noise parameters
    if radar_hparams.algo == vars.imm_name
        radar_hparams.Q{1} = eye(3) * 30^2;
        radar_hparams.Q{2} = eye(3) * 80^2;
        radar_hparams.Q{3} = eye(3) * 200^2;
        radar_hparams.n_mode = length(radar_hparams.Q);

        val = 0.9; % Probability that we stay in the same mode
        radar_hparams.PI = ones(radar_hparams.n_mode) .* (1-val) ./ (radar_hparams.n_mode - 1); % Rest of the transitions are uniform
        for i = 1:radar_hparams.n_mode
            radar_hparams.PI(i, i) = val;
        end
    else
        radar_hparams.Q = eye(3) * 100^2;
    end

    %% Sensor noise parameters
    radar_hparams.R = diag([50; 0.01; 0.01])^2;
    % radar_hparams.R = diag([50; 0.05; 0.05])^2;

    %% Covariance initalization parameters
    radar_hparams.v_max = 600;
    radar_hparams.kappa = 2;

    %% Azimuth modding
    radar_hparams.modding_param_index = 2;

    %% A and B matrices
    radar_hparams.A = const_vel3;
    radar_hparams.B = B3;
    radar_hparams.init_vel = [0, 0, 0]';

    %% h and inverse h 
    radar_hparams.inv_h = @(y) sph2Cart_(y);
    radar_hparams.h     = @(x) cart2Sph_(x);
    
    %% Parameter selection for particle filter
    if radar_hparams.algo == vars.pf_name
        radar_hparams.n_particles = 5000;
        radar_hparams.resample_method = "Multinomial";
        radar_hparams.need_resample_method = "MaxWeight";
        radar_hparams.resample_parameter = 5 / radar_hparams.n_particles;
    end
end
%%
sync_measurements = false;       % Whether all sensors measure at exactly the same time

%% Os details
os_color = [0, 0, 0];
os_init_pos = [0, 0, 0]';
% os_init_pos = [6000, 6000, 6000]';

%% Target details
ref_pos= [-200000 200000; -200000 200000; -10000 10000];
ref_vel = [-300 300; -300 300; -10 10];
%%
process_bias = [0, 0, 0];
process_bias_os = [0, 0, 0];

Q_target = eye(3) * 100^2;
Q_os     = eye(3) * 100^2;
% Q_os     = eye(3) * 0^2;

%% Initial positions and velocity
tgt_init_pos =  2 * (rand(3, n_target) - 0.5 * ones(3,1)) .* (diff(ref_pos, 1, 2) - ref_pos(:, 2));
% tgt_init_pos = [1.4 -0.55 -0.9 1.9 0.2;
%                 -1.2 -0.9 0.1 1.9 2.4;
%                 -0.01 0.006 -0.0777 0.06 -0.042]*1e5;

tgt_init_vel =  2 * (rand(3, n_target) - 0.5 * ones(3,1)) .* (diff(ref_vel, 1, 2) - ref_vel(:, 2));
os_init_vel  =  2 * (rand(3, 1) - 0.5 * ones(3,1)) .* (diff(ref_vel, 1, 2) - ref_vel(:, 2));
% os_init_vel  =  [0 0 0]';

%% Cell initalization
tgt_state = cell(1, n_target);
meas = struct();
os = struct();

tgt_meas_times = cell(n_target, n_sensor);
tgt_state_meas_time = cell(n_target, n_sensor);
os_state_meas_time = cell(n_target, n_sensor);
tgt_state_rel_meas_time = cell(n_target, n_sensor);

os.state = zeros(6, length(sim_t_arr));

if include_radar
    for m = 1:n_monte
        meas(m).radar.sph  = cell(1, n_target);
        meas(m).radar.cart = cell(1, n_target);
    end
end

tgt_state_times = cell(1, n_target);
fa_times = cell(1, n_sensor);

%% True state propagation
for i = 1:n_target
    tgt_state{i} = zeros(6, length(sim_t_arr));
    tgt_state{i}(:, 1) = [tgt_init_pos(:, i); tgt_init_vel(:, i)];
end

for i = 1:n_target
    for k = 2:length(sim_t_arr)
        tgt_state{i}(:, k) = global_motion_model * tgt_state{i}(:, k - 1) + B3(sim_dt) * mvnrnd_NT(process_bias, Q_target)';
    end
end

for i = 1:n_target
   tgt_state_times{i} = sim_t_arr;
end

os.state(:, 1) = [os_init_pos; os_init_vel];
for k = 2:length(sim_t_arr)
    os.state(:, k) = global_motion_model * os.state(:, k - 1) + B3(sim_dt) * mvnrnd_NT(process_bias_os, Q_os)';
end

rng(time_seed)
%% Measurement
for i = 1:n_target
    for j = 1:n_sensor
        sensor_name = sensor_list(j);
        if sensor_name == vars.radar_name
            t_mean = radar_hparams.t_mean;
            t_var  = radar_hparams.t_var;
        else
            sensorNameError(sensor_name)
        end
        if ~sync_measurements || j < 2
            tgt_meas_times{i, j} = getMeasurementTimes(sim_t_arr(1), sim_t_arr(end), t_mean, t_var);
        else
            tgt_meas_times{i, j} = tgt_meas_times{i, j-1};
        end        
    end
end

if is_seeded
    rng(other_seed)
else
    rng('shuffle')
end

for m = 1:n_monte
    for i = 1:n_target
        for j = 1:n_sensor
            sensor_name = sensor_list(j);
            if sensor_name == vars.radar_name
                meas(m).radar.time{i} = tgt_meas_times{i, j};
            else
                sensorNameEror(sensor_name)
            end
        end
    end
end

for i = 1:n_target
    for j = 1:n_sensor
        tgt_state_meas_time{i, j} = interp1(sim_t_arr, tgt_state{i}', tgt_meas_times{i, j})';
    end
end

for i = 1:n_target
    for j = 1:n_sensor
        os_state_meas_time{i, j} = interp1(sim_t_arr, os.state', tgt_meas_times{i, j})';
    end
end

for i =1:n_target
    for j = 1:n_sensor
        tgt_state_rel_meas_time{i, j} = tgt_state_meas_time{i, j} - os_state_meas_time{i, j};
    end
end

noises = cell(n_target, n_sensor);
for m = 1:n_monte
    for i = 1:n_target
        for j = 1:n_sensor
            sensor_name = sensor_list(j);
            if sensor_name == vars.radar_name
                meas_true_radar          = cart2Sph_(tgt_state_rel_meas_time{i, j});
                noises{i, j} = mvnrnd_NT(radar_hparams.bias ...
                                                , radar_hparams.R, length(meas_true_radar))';
                meas(m).radar.sph{i} = meas_true_radar + noises{i, j};
                meas(m).radar.sph{i}(2:3, :) = azElModder(meas(m).radar.sph{i}(2:3, :));
                meas(m).radar.cart{i}       = sph2Cart_(meas(m).radar.sph{i});
            else
                sensorNameError(sensor_name)
            end
        end
    end
end

for m = 1:n_monte
    for j = 1:n_sensor
        sensor_name = sensor_list(j);
        if sensor_name == vars.radar_name
            meas(m).radar.sph_all = cell2Mat_(meas(m).radar.sph);
        else
            sensorNameError(sensor_name)
        end
    end
end

for m = 1:n_monte
    for j = 1:n_sensor
        sensor_name = sensor_list(j);
        if sensor_name == vars.radar_name
            meas(m).radar.time{i} = tgt_meas_times{i, j};
            [meas(m).radar.all_times, sort_idx] = sort(cell2Mat_(meas(m).radar.time));
            meas(m).radar.sph_all = meas(m).radar.sph_all(:, sort_idx);
            id_sizes = zeros(1, n_target);
            for i = 1:n_target
                id_sizes(i) = length(meas(m).radar.time{i});
            end
            meas(m).radar.id = [];
            for i = 1:n_target
                meas(m).radar.id = [meas(m).radar.id repmat(i, 1, id_sizes(i))];
            end
            meas(m).radar.id = meas(m).radar.id(sort_idx);
        else
            sensorNameEror(sensor_name)
        end
    end
end

os.state_meas_time = struct();
for m = 1:n_monte
    for j = 1:n_sensor
        sensor_name = sensor_list(j);
        if sensor_name == vars.radar_name
            os.state_meas_time(j).sensor_name =  sensor_name;
            os.state_meas_time(j).state = interp1(sim_t_arr, os.state', meas(m).radar.all_times)';
            os.state_meas_time(j).time = meas(m).radar.all_times;
        else
            sensorNameError(sensor_name)
        end
    end
end

%% False alarm
measurement_id_fa = [];
os.state_meas_time_w_fa = os.state_meas_time;
for m = 1:n_monte
    for j = 1:n_sensor
        sensor_name = sensor_list(j);
        if sensor_name == vars.radar_name
            if radar_hparams.lambda == 0
                nfa = 0;
            else
                nfa = poissrnd_NT(radar_hparams.lambda);
            end
            fa_times{j} = calculateFATimes(sim_t_arr(1), sim_t_arr(end), nfa, radar_hparams.lambda);
            fa_indexes = bsLess(meas(m).radar.all_times, fa_times{j});
            %% TODO: Reduce this fa state to FOV
            fa_state =  [2 * (rand(3, length(fa_indexes)) - 0.5 * ones(3,1)) .* (diff(ref_pos, 1, 2) - ref_pos(:, 2));
                 2 * (rand(3, length(fa_indexes)) - 0.5 * ones(3,1)) .* (diff(ref_vel, 1, 2) - ref_vel(:, 2))];
            fa_state_sph = cart2Sph_(fa_state);
            meas(m).radar.sph_all_w_fa = insertFA(meas(m).radar.sph_all, fa_state_sph, fa_indexes);
            meas(m).radar.all_times_w_fa = insertFA(meas(m).radar.all_times, fa_times{j}, fa_indexes);
            meas(m).radar.id_w_fa = insertFA(meas(m).radar.id, -1*ones(1, length(fa_indexes)), fa_indexes);
            meas(m).radar.fa_indexes = fa_indexes + (1:length(fa_indexes));
            if nfa ~= 0
                os_state_fa_times = interp1(sim_t_arr, os.state', fa_times{j})';
                os.state_meas_time_w_fa(j).state = insertFA(os.state_meas_time_w_fa(j).state, os_state_fa_times, fa_indexes);
                os.state_meas_time_w_fa(j).time = insertFA(os.state_meas_time_w_fa(j).time, fa_times{j}, fa_indexes);
            end
            clear os_state_fa_times
        else
            sensorNameError(sensor_name)
        end
    end
end

%% Batch tracking (Radar)
if include_radar
    radar_logs = cell(1, n_monte);
    for m = 1:n_monte
        if radar_hparams.algo == vars.ekf_name || radar_hparams.algo == vars.ukf_name
            [radar_logs{m}] = radarTracker(meas(m).radar, radar_hparams, vars);
        elseif radar_hparams.algo == vars.pf_name
            [radar_logs{m}] = radarTrackerPF(meas(m).radar, radar_hparams, vars);
        elseif radar_hparams.algo == vars.imm_name
            [radar_logs{m}] = radarTrackerIMM(meas(m).radar, radar_hparams, vars);            
        else
            error("This algorithm name for Radar is not handled")
        end
    end
end

max_num_radar = 0;
for m = 1:n_monte
    if ~isempty(radar_logs{m})
        for i = 1:length(radar_logs{m}.confirmed_track_log)
            max_num_radar = max(max_num_radar, length(radar_logs{m}.confirmed_track_log{i}));
        end
    end
end


%%
% save('data.mat')
% load('data.mat')
%%


colors = chooseFromColorMap([], n_target);
figure(); hold on; grid on; view(3)
plot3(os.state(1, :), os.state(2, :), os.state(3, :), 'Color', os_color, 'LineWidth',2);
for i = 1:n_target
    plot3(tgt_state{i}(1, :)...
    , tgt_state{i}(2, :) ...
    , tgt_state{i}(3, :)...
    ,'Color', colors(i, :), 'LineWidth', 4) 
end

max_tgt_index_mc = 0;
for i = 1:length(radar_logs{index_mc}.confirmed_track_log)
    for j = 1:length(radar_logs{index_mc}.confirmed_track_log{i})
        max_tgt_index_mc = max(max_tgt_index_mc, radar_logs{index_mc}.confirmed_track_log{i}(j).track_id);
    end
end

len_track_index_mc = length(radar_logs{index_mc}.confirmed_track_log);
state_cell_radar = cell(1, max_tgt_index_mc);
time_cell_radar = cell(1, max_tgt_index_mc);
for i = 1:max_tgt_index_mc
    state_cell_radar{i} = [];
    time_cell_radar{i} = [];
end

for i = 1:len_track_index_mc
    for j = 1:length(radar_logs{index_mc}.confirmed_track_log{i})
        id = radar_logs{index_mc}.confirmed_track_log{i}(j).track_id;
        time_index = max([radar_logs{index_mc}.confirmed_track_log{i}.t]);
        state_cell_radar{id}(:, end+1) = radar_logs{index_mc}.confirmed_track_log{i}(j).x;
        time_cell_radar{id}(end+1) = time_index;
    end
end

times_radar = meas(index_mc).radar.all_times_w_fa;
index_os = zeros(1, length(times_radar));
for j = 1:length(times_radar)
    val = find(os.state_meas_time_w_fa.time == times_radar(j));
    assert(~isempty(val));
    val = val(1);
    index_os(j) = val;
end
%% Plot radar measurements
cart_pos_radar = sph2Cart_(meas(index_mc).radar.sph_all_w_fa);

x_vals_radar = cart_pos_radar(1, :);
y_vals_radar = cart_pos_radar(2, :);
z_vals_radar = cart_pos_radar(3, :);

x_os = os.state_meas_time_w_fa.state(1, index_os);
y_os = os.state_meas_time_w_fa.state(2, index_os);
z_os = os.state_meas_time_w_fa.state(3, index_os);

scatter3(x_vals_radar + x_os, ...
         y_vals_radar + y_os, ...
         z_vals_radar + z_os, ...
         5, ...
         'MarkerEdgeColor', [0 0 0]', ...
         'MarkerFaceColor', [0 0 0]');

%% Plot estimations
for i = 1:max_tgt_index_mc
    x_vals_radar_est = state_cell_radar{i}(1, :);
    y_vals_radar_est = state_cell_radar{i}(2, :);
    z_vals_radar_est = state_cell_radar{i}(3, :);

    times_radar_est = time_cell_radar{i};

    index_os = zeros(1, length(times_radar_est));
    for j = 1:length(times_radar_est)
        val = find(os.state_meas_time_w_fa.time == times_radar_est(j));
        assert(~isempty(val));
        val = val(1);
        index_os(j) = val;
    end

    x_os = os.state_meas_time_w_fa.state(1, index_os);
    y_os = os.state_meas_time_w_fa.state(2, index_os);
    z_os = os.state_meas_time_w_fa.state(3, index_os);

    scatter3(x_vals_radar_est + x_os, ...
             y_vals_radar_est + y_os, ...
             z_vals_radar_est + z_os, ...
             5, 'filled', ...
             'MarkerEdgeColor', colors(i, :), ...
             'MarkerFaceColor', colors(i, :), ...
             'Marker','*');
end

title("Monte Carlo Run " + num2str(index_mc))
disp("Plotted mc run " + num2str(index_mc))

