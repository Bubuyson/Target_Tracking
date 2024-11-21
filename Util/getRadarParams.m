function radar_hparams = getRadarParams(algo_name, vars)

    const_vel3 = @(dt) [eye(3) eye(3) * dt; zeros(3) eye(3)];
    B3 = @(dt) [eye(3) * dt^2/2; eye(3)* dt];

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
    radar_hparams.algo = algo_name;

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
    % radar_hparams.R = diag([50; 0.001; 0.001])^2;
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
    if algo_name == vars.pf_name
        radar_hparams.n_particles = 5000;
        radar_hparams.resample_method = "Multinomial";
        radar_hparams.need_resample_method = "MaxWeight";
        radar_hparams.resample_parameter = 5 / radar_hparams.n_particles;
    end


end