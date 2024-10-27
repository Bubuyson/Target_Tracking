clc; clear; close all
%%
% rng(6)
%%
T=1;
%%
I2=eye(2);
Z2=zeros(2);
A = [I2 T*I2 ; Z2 I2];
B =[T^2/2*I2 ; T*I2];
%%
x0_M = [1000,1000,100,100];
P0 = diag([100,100,10,10])^2;
x = mvnrnd(x0_M,P0,1)';
x0 = x';
n=100;
n_particles = 2000;
%%
need_resample_method = "MaxWeight"; max_w = 5 / n_particles;
% need_resample_method = "EffectiveSize"; n_eff = (2/3) * n_particles;
resample_method_name = "Multinomial";
% resample_method_name = "Stratified";
% resample_method_name = "Systematic";
% resample_method_name = "Residual";
%%
Q=eye(2)*100;
R=diag([10^2,(25*pi/180)^2]); 
%%

for i=2:n
   x(:,i) = A*x(:,i-1)+B*mvnrnd(zeros(1,2),Q,1)';
end

h=@(x) [sqrt(x(1,:).^2 + x(2,:).^2); atan2(x(2,:),x(1,:))]; 
y=h(x)+mvnrnd(zeros(1,2),R,n)';

f =@(T) [I2 T*I2 ; Z2 I2];
B =@(T) [T^2/2*I2 ; T*I2];

mean_particles_log = cell(1, n);
particles_log = cell(1, n);

particles_0 = mvnrnd(x(:, 1), P0, n_particles)';
weights = ones(1, n_particles) ./ n_particles;
mean_particles_log{1} = particles_0 * weights';
particles_log{1} = particles_0;
prev_particles = particles_0;
for i = 2:n
    %% State transition
    new_particles = f(T) * prev_particles + B(T) * mvnrnd([0 0]', Q, n_particles)';
    %% Calculate likelihood
    y_ = h(new_particles);
    meas = y(:, i);
    likelihoods = mvnpdf(meas', y_', R)';
    weights = weights .* likelihoods;
    weights = weights ./ sum(weights);
    %% Resample if necessary
    indices = 1:n_particles;
    if need_resample_method == "MaxWeight"
        if max(weights) > max_w
            indices = do_resample(weights, resample_method_name);
            weights = ones(1, n_particles) ./ n_particles;
        end
    elseif need_resample_method == "EffectiveSize"
        if 1/sum(weights.^2) < n_eff
            indices = do_resample(weights ,resample_method_name);
            weights = ones(1, n_particles) ./ n_particles;
        end        
    else
        error("Name error for need resampling method")
    end
    new_particles = new_particles(:, indices);  
    %% Calculate mean particle (and covariance if you want detailed analysis)
    mean_particle = new_particles * weights';
    mean_particles_log{i} = mean_particle;
    particles_log{i} = new_particles;
    prev_particles = new_particles;
end

figure(); hold on; grid on;
for i = 1:n
    scatter(particles_log{i}(1, :), particles_log{i}(2, :))
    scatter(mean_particles_log{i}(1), mean_particles_log{i}(2), 'r*')
end
plot(x(1, :), x(2, :), 'k', 'linewidth', 2)
xlabel('x'); ylabel('y'); title('Trajectories')

function indices = do_resample(weights, resample_method_name)
    if resample_method_name == "Systematic"
        indices = systematic_resampling(weights);
    end
    if resample_method_name == "Residual"
        indices = residual_resampling(weights);
    end
    if resample_method_name == "Stratified"
        indices = stratified_resampling(weights);
    end
    if resample_method_name == "Multinomial"
        indices = multinomial_resampling(weights);
    end

    function indices = systematic_resampling(weights)

    end

    function indices = residual_resampling(weights)

    end

    function indices = stratified_resampling(weights)

    end

    function indices = multinomial_resampling(weights)
        n = length(weights);
        r = rand(1, n);
        cdf = cumsum(weights);
        indices = zeros(1, n);
        for i = 1:n
            indices(i) = find(r(i) < cdf, 1);
        end
    end
end




