classdef ParticleFilter

    properties (Access = public)
        particles
        mean_particle
        mean_particle_log
        num_particle
        particle_log
        weights
        f
        B
        Q
        h
        R
        last_dt
        num_state
        num_meas
%         cycle = 0
        time_log = [0];
        resample_method = "Multinomial"         % methods = ["Multinomial", "Residual", "Systematic", "Stratified"];
        resample_function = "MaxWeight"         % functions = ["MaxWeights", "EffectiveSize", "Every"];
        resample_function_parameter = 5;
        apply_MH = false;
        apply_auxiliary = false;
        log = true;
    end

    properties (Access = protected)
        every_counter = 0;
        auxiliary_weights;
    end
    
    methods (Access = public)
        function obj = ParticleFilter(num_particle, f, B, Q, R, h)
            if nargin < 3
                obj = obj.setDefaults();
            else
                assert(isdiag(Q), "Input a diagonal positive definite process noise covariance")
                assert(isdiag(R), "Input a diagonal positive definite measurement noise covariance")
                obj.f = f;
                obj.B = B;
                obj.Q = Q;
                obj.h = h;
                obj.R = R;
            end
            obj.num_state = length(obj.f(0));   % Since Q may be defined as acceleration variance rather than velocity and position variance
            obj.num_meas = length(obj.R);
            obj.particles = NaN*ones(obj.num_state, num_particle);
            obj.weights = ones(1, num_particle) / num_particle;
            obj.num_particle = num_particle;
        end

        function obj = setDefaults(obj)
            I2=eye(2);
            Z2=zeros(2);
            obj.f = @(dt) [I2 dt*I2 ; Z2 I2];
            obj.B = @(dt) [dt^2/2*I2 ; dt*I2];
            obj.Q = eye(2)*100;
            obj.R = diag([10^2,(25*pi/180)^2]); 
            obj.h = @(x) [sqrt(x(1,:).^2 + x(2,:).^2); atan2(x(2,:),x(1,:))];
        end

        function obj = setf(obj, f)
            obj.f = f;
        end

        function obj = setQ(obj, Q)
            obj.Q = Q;
        end

        function obj = setB(obj, B)
            obj.B = B;
        end

        function obj = seth(obj, h)
            obj.h = h;
        end

        function obj = setR(obj, R)
            obj.R = R;
        end

        function obj = setResample(obj, resample_method, varargin)
            methods = ["Multinomial", "Residual", "Systematic", "Stratified"];
            functions = ["MaxWeight", "EffectiveSize", "Every"];
            if nargin == 2 && ~isempty(obj.isIn(methods, resample_method))
                obj.resample_method = resample_method;

            elseif nargin == 3 && ~isempty(obj.isIn(methods, resample_method)) ...
                && ~isempty(obj.isIn(functions, varargin{1}))
                obj.resample_method = resample_method;
                obj.resample_function = varargin{1};

            elseif nargin == 4 && ~isempty(obj.isIn(methods, resample_method)) ...
                    && ~isempty(obj.isIn(functions, varargin{1})) ...
                    && isa(varargin{2}, 'double')
                obj.resample_method = resample_method;
                obj.resample_function = varargin{1};
                obj.resample_function_parameter = varargin{2};
            else
                error("Input a valid resample method (and resample function and its ratio if you want)")
            end
        end

        function obj = setApplyMH(obj, flag)
            assert(isa(flag, 'logical'), "Enter true or false for MH")
            obj.apply_MH = flag;
        end

        function obj = setApplyAuxiliary(obj, flag)
            assert(isa(flag, 'logical'), "Enter true or false for Auxiliary Particle Filter")
            obj.apply_auxiliary = flag;
        end

        function obj = setLog(obj, flag)
            assert(isa(flag, 'logical'), "Enter true or false for logging")
            obj.log = flag;
        end

        function obj = initialize(obj, x0, P0)
            obj.particles = obj.mvnrnd_NT(x0, P0, obj.num_particle)';
        end

        function obj = update(obj, y, dt)
            
            if obj.apply_auxiliary
                obj = obj.auxiliaryUpdate(y, dt);
            elseif ~obj.apply_auxiliary
                obj = obj.standardUpdate(y, dt);
            end
            
        end

        function obj = standardUpdate(obj, y, dt)
            obj.last_dt = dt;
            obj = obj.noisyStateTransition();
            obj.weights = obj.weights .* obj.likelihoodCalculator(y)';
            obj.weights = obj.weights / sum(obj.weights);

            if obj.needResampling()
                indices = obj.resample();
                obj.weights = ones(1, obj.num_particle) / obj.num_particle;
                if obj.apply_MH
                    candidate_particles = obj.particles(:, indices);
                    obj = obj.metropolisHasting(candidate_particles, y);
                else
                    obj.particles = obj.particles(:, indices);
                end
            end
            obj.mean_particle = obj.particles * obj.weights';
            if obj.log
                obj.particle_log(:, :, end+1) = obj.particles;
                obj.mean_particle_log(:, end+1) = obj.mean_particle;
                obj.time_log(end+1) = obj.time_log(end) + dt;
            end
            obj.last_dt = dt;
        end

        function obj = auxiliaryUpdate(obj, y, dt)
            obj.last_dt = dt;
            obj = obj.noisyStateTransition();
            obj.auxiliary_weights = obj.weights .* obj.likelihoodCalculator(y)';
            obj.auxiliary_weights = obj.auxiliary_weights / sum(obj.auxiliary_weights);

           
            obj.mean_particle = obj.particles * obj.weights';
            if obj.needResampling()
                indices = obj.resample();
                obj.weights = ones(1, obj.num_particle) / obj.num_particle;
                if obj.apply_MH
                    candidate_particles = obj.particles(:, indices);
                    obj = obj.metropolisHasting(candidate_particles, y);
                else
                    obj.particles = obj.particles(:, indices);
                end
            end
            
            
            if obj.log
                obj.particle_log(:, :, end+1) = obj.particles;
                obj.mean_particle_log(:, end+1) = obj.mean_particle;
                obj.time_log(end+1) = obj.time_log(end) + dt;
            end
            obj.last_dt = dt;
        end

        function plotTrajectory(obj, x_true, dims)
            if ~obj.log
                error("Enable logging for plotting trajectory")
            end
            Ndims = length(dims);
            figure(); hold on; grid on;
            if Ndims == 2
                for i = 1:size(obj.particle_log, 3)
                    scatter(obj.particle_log(dims(1), :, i), obj.particle_log(dims(2), :, i))
                end
                plot(x_true(dims(1), :), x_true(dims(2), :), 'k', 'linewidth', 2)
                scatter(obj.mean_particle_log(dims(1), :), obj.mean_particle_log(dims(2), :), 'r*')
                xlabel('x'); ylabel('y'); title('Trajectories')
            elseif Ndims == 3
                for i = 1:size(obj.particle_log, 3)
                    scatter3(obj.particle_log(dims(1), :, i), obj.particle_log(dims(2), :, i)...
                        , obj.particle_log(dims(3), :, i))
                end
                plot3(x_true(dims(1), :), x_true(dims(2), :), x_true(dims(3), :), 'k', 'linewidth', 2)
                scatter3(obj.mean_particle_log(dims(1), :), obj.mean_particle_log(dims(2), :)...
                    , obj.mean_particle_log(dims(3), :), 'r*')
                xlabel('x'); ylabel('y'); zlabel(z); title('Trajectories')
            else
                error("Input correct number of dimensions")
            end
        end

        function calculateRMSE(obj, x_true, dims, plot_option)
            if ~obj.log
                error("Enable logging for calculating and plotting rmse")
            end
            if nargin == 3
                plot_option = false;
            elseif nargin == 4 
            else
                error("Input true trajectory and dimensions (and plot flag if desired)")
            end
            Ndims = size(dims);
            calculate_velocity = Ndims(1) == 2;
            pos_rmse = sqrt(sum((x_true(dims(1, :), :) - obj.mean_particle_log(dims(1, :), :)).^2));
            mean_rms_pos = mean(pos_rmse);
            disp("The average RMS position error of the particle filter is " + int2str(mean_rms_pos))
            if calculate_velocity
                vel_rmse = sqrt(sum((x_true(dims(2, :), :) - obj.mean_particle_log(dims(2, :), :)).^2));
                mean_rms_vel = mean(vel_rmse);
                disp("The average RMS velocity error of the particle filter is " + int2str(mean_rms_vel))
            end
            
            if plot_option
                figure(); hold on; grid on;
                subplot(1 + calculate_velocity, 1, 1)
                plot(obj.time_log(2:end), pos_rmse); grid on;
                xlabel("Time")
                ylabel("Position RMSE")
                title("Position RMSE with respect to time")
                if calculate_velocity
                    subplot(1 + calculate_velocity, 1, 2)
                    plot(obj.time_log(2:end), vel_rmse); grid on;
                    xlabel("Time")
                    ylabel("Velocity RMSE")
                    title("Velocity RMSE with respect to time")
                end
            end
        end

    end

    methods (Access = protected)
        function obj = metropolisHasting(obj, x_candidate, y)
            % TODO: Work on parallelization (it should be possible)
            N = obj.num_particle;
            for i = 1:N
                x_candidate_i = x_candidate(:, i);
                x_current_i = obj.particles(:, i);
                acceptance_probability = min(1, obj.likelihoodCalculatorExternal(y, x_candidate_i )/...
                                            obj.likelihoodCalculatorExternal(y, x_current_i)); % candidate over current
                
                if rand() <= acceptance_probability     % Choose to accept or reject the update
                    obj.particles(:, i) = x_candidate(:, i);
                end
            end
        end


        function obj = noisyStateTransition(obj)
            N = obj.num_particle;
            obj.particles = obj.f(obj.last_dt) * obj.particles...
                + obj.B(obj.last_dt) * mvnrnd_NT(zeros(1, 2), obj.Q, N)';
        end

        function x = noisyStateTransitionCalculateOne(obj, index)
            N = obj.num_particle;
            element = obj.particles(:, index);
            x = obj.f(obj.last_dt) * element + obj.B * mvnrnd_NT(zeros(1, 2), obj.Q, N)';

        end

        function likelihood = likelihoodCalculator(obj, y)
            sigma = obj.R;
            y_ = obj.h(obj.particles);
            [n_state, ~] = size(y);
            [n_state_2, n_particle] = size(y_);
            assert(n_state == n_state_2, "Error in x and mu dimensions")
            likelihood = zeros(n_particle, 1);
            for i = 1: n_particle
                likelihood(i) = 1/sqrt((2*pi)^length(y)*det(sigma))*exp(-1/2*(y-y_(:, i))'/sigma*(y-y_(:, i)));
            end
        end

        function likelihood = likelihoodCalculatorExternal(obj, y, x_)
            sigma = obj.R;
            y_ = obj.h(x_);
            [n_state, ~] = size(y);
            [n_state_2, n_particle] = size(y_);
            assert(n_state == n_state_2, "Error in x and mu dimensions")
            likelihood = zeros(n_particle, 1);
            for i = 1: n_particle
                likelihood(i) = 1/sqrt((2*pi)^length(y)*det(sigma))*exp(-1/2*(y-y_(:, i))'/sigma*(y-y_(:, i)));
            end
        end

        function [flag, obj]= needResampling(obj)
            switch obj.resample_function
                case "MaxWeight"
                    flag = obj.isOverMaxWeight();
                case "EffectiveSize"
                    flag = obj.isOverEffectiveSize();
                case "Every"
                    [flag, obj] = obj.isEvery();
            end
        end

        function flag = isOverEffectiveSize(obj)
            if obj.apply_auxiliary
                N_eff = calculateEffectiveSize(obj.auxiliary_weights);
            else
                N_eff = calculateEffectiveSize(obj.weights);
            end
            flag = N_eff < obj.num_particle / obj.resample_function_parameter;
            function N_eff = calculateEffectiveSize(weights)
                N_eff = 1/sum(weights.^2);
            end
        end

        function flag = isOverMaxWeight(obj)
            if obj.apply_auxiliary
                maxval = max(obj.auxiliary_weights);
            else
                maxval = max(obj.weights);
            end
            flag = maxval > obj.resample_function_parameter / obj.num_particle;
        end

        function [flag, obj] = isEvery(obj)
            flag = false;
            obj.every_counter = obj.every_counter + 1;
            if mod(obj.every_counter, obj.resample_function_parameter)
                flag = true;
                obj.every_counter = 0;
            end
        end

        function chosen_indices = resample(obj)
            % it is said that residual and stratified resampling outperform
            % multinomial  resampling (for sure) and systematic resampling 
            % (sometimes) [Elfring, p.16]
        
            %methods = ["Multinomial", "Residual", "Systematic", "Stratified"];
            if obj.apply_auxiliary
                weights_ = obj.auxiliary_weights;
            else
                weights_ = obj.weights;
            end
            switch obj.resample_method
                case "Systematic"
                    chosen_indices = obj.systematicResample(weights_);
                case "Multinomial"
                    chosen_indices = obj.multinomialResample(weights_);
                case "Residual"
                    chosen_indices = obj.residualResample(weights_);
                case "Stratified"
                    chosen_indices = obj.stratifiedResample(weights_);
            end
    
        end

        function chosen_indices = systematicResample(obj, weights)
            % The same as stratified resample, where the only difference is
            % the random number is only generated once
            N = obj.num_particle;
            rand_vals = rand(1);        % THIS ONE IS GENERATED ONCE
            positions = (linspace(0, N-1, N) + rand_vals) / N;
            cum_sum = cumsum(weights);
            cum_sum(end) = 1;
            i = 1;  % tracks the current selection point
            j = 1;  % tracks the current particle
            chosen_indices = zeros(1, N);
            while i <= N
                if positions(i) < cum_sum(j)
                    chosen_indices(i) = int32(j);
                    i = i+1;
                else
                    j = j+1;
                end
            end
        end
        
        function chosen_indices = multinomialResample(obj, weights)
            % Generates N uniform random numbers and puts them to the y axis of
            % cdf of weights and takes corresponding particle in x axis
            N = obj.num_particle;
            cum_sum = cumsum(weights);
            cum_sum(end) = 1;
            chosen_indices = zeros(1, N);
            for i = 1:N
                rand_val = rand(1);
                chosen_indices(1, i) = obj.binarySearchLess(cum_sum, rand_val);            % binary search over the cdf
%                 chosen_indices(1, i) = find(rand_val < cum_sum, 1);
            end
        end
    
        function chosen_indices = residualResample(obj, weights)
            % Firstly generates some particles using the floor function,
            % and the rest of the particles are chosen via multinomial
            % resampling
            N = obj.num_particle;
            floors = floor(N * weights);
            k = 1;
            chosen_indices = zeros(1, N);
            for i = 1:N
                for j = 1:floors(i)
                    chosen_indices(k) = i;
                    k = k+1;
                end
            end
            residual = weights - floors;
            residual = residual / sum(residual);
            cum_sum = cumsum(residual);
            cum_sum(1) = 0;
            cum_sum(end) = 1;
            for i = k:N
                rand_val = rand();
                chosen_indices(i) = obj.binarySearchLess(cum_sum, rand_val);
%                 chosen_indices(i) = find(cum_sum < rand_val, 1);
            end
        end
    
        function chosen_indices = stratifiedResample(obj, weights)
            % Generates ordered random numbers and does multinomial
            % resample using these numbers
            N = obj.num_particle;
            positions = (rand(1, N) + linspace(0, N-1, N)) / N;
            cum_sum = cumsum(weights);
            cum_sum(end) = 1;
            i = 1;
            j = 1;
            chosen_indices = zeros(1, N);
            while i <= N
                if positions(i) < cum_sum(j)
                    chosen_indices(i) = j;
                    i = i + 1;
                else
                    j = j + 1;
                end
            end
        end
    end

    % Utility
    methods (Access = private)

        function output = isIn(~, array, val)
            output = find(array == val);
        end

        function index = binarySearchLess(~, sortedArray, target)
            % Returns the index of the first element that is less than the searched
            % value, assuming target can not be more than the last element of the
            % array
            l = int32(1); r = int32(length(sortedArray));
            assert(target <= sortedArray(end), 'Choose target less than the maximum (end) value of the array')
            while l<=r
                m = (l + r) / 2;
                if(sortedArray(m) == target)
                    index = m;
                    break
                end
                if(sortedArray(m) > target)
                    r = m-1;
                    index = m;
                else
                    l = m+1;
                end
            end

        end

        function out = mvnrnd_NT(~, mu, sigma, num)
            sq_sig = sqrt(sigma);
            out = repmat(mu,num,1) + randn(num, size(mu, 2))*sq_sig;
        end

    end

end