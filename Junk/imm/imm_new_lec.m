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
%%
Q={eye(2)*100; eye(2)*1000; eye(2)*10};
R=diag([10^2,(25*pi/180)^2]); 

%%
for i=2:n
   x(:,i) = A*x(:,i-1)+B*mvnrnd(zeros(1,2),Q{2},1)';
end

h=@(x) [sqrt(x(1,:).^2 + x(2,:).^2); atan2(x(2,:),x(1,:))]; 
y=h(x)+mvnrnd(zeros(1,2),R,n)';

f =@(T) [I2 T*I2 ; Z2 I2];
B =@(T) [T^2/2*I2 ; T*I2];


x_log       =   cell(1, n);
x_mean_log  =   cell(1, n);
P_log       =   cell(1, n);
P_mean_log  =   cell(1, n);
k_log       =   cell(1, n);
%%
PI = [0.90 0.05 0.05
      0.05 0.90 0.05
      0.05 0.05 0.90];

%%
imm = struct();
imm(1).state.x = x0';
imm(1).state.P = P0;
imm(1).state.k = 1/3;
imm(1).f = A;
imm(1).B = B;
imm(1).Q = Q{1};

imm(2).state.x = x0';
imm(2).state.P = P0;
imm(2).state.k = 1/3;
imm(2).f = A;
imm(2).B = B;
imm(2).Q = Q{2};

imm(3).state.x = x0';
imm(3).state.P = P0;
imm(3).state.k = 1/3;
imm(3).f = A;
imm(3).B = B;
imm(3).Q = Q{3};

imm(1).x_out(:, 1) = x0';
imm(1).P_out(:, :, 1) = P0;

n_mode = length(imm);

for i = 1:n_mode
    x_log{1}(:, i)    = imm(i).state.x;
    P_log{1}(:, :, i) = imm(i).state.P;
    k_log{1}(i)       = imm(i).state.k;
end

x_mean_log{1} = x_log{1} * k_log{1}';
for i = 2:n

    %% Predict weights
    mu = [imm(1).state.k imm(2).state.k imm(3).state.k]';
    mu_ = (PI)^T * mu;
    MU = (PI)^T .* (repmat(mu',n_mode,1)) ./ repmat(mu_,1,n_mode); 
    %% Merging
    for j=1:n_mode
        x_IMM(:, 1)    = imm(1).state.x;
        x_IMM(:, 2)    = imm(2).state.x;
        x_IMM(:, 3)    = imm(3).state.x;

        P_IMM(:, :, 1) = imm(1).state.P;
        P_IMM(:, :, 2) = imm(2).state.P;
        P_IMM(:, :, 3) = imm(3).state.P;

        [xhat(:,j), Phat(:, :, j)] = mergeIMM(MU(j,:), x_IMM, P_IMM);
    end

    %% State transition
    for j = 1:n_mode
        x_(:, j) = f(T) * imm(j).state.x;
        P_(:, :, j) = f(T) * imm(j).state.P * f(T)' + B(T) * Q{j} * B(T)';
    end

    %% Calculate likelihood
    y_ = h(x_);
    meas = y(:, i);
    for j=1:n_mode
        H = getNumericalJacobian(h, x_(:,j));
        S = (H * P_(:, :, j) * H' + R);
        K = P_(:, :, j) * H' / S;
        
        imm(j).state.x = x_(:, j) + K * (meas - y_(:, j));
        imm(j).state.P = P_(:, :, j) - K * S * K';
        imm(j).state.P = (imm(j).state.P + imm(j).state.P') ./ 2;

        likelihood(j) = mvnpdf(meas, y_(:,j), S);

        
    end
    k = (likelihood .* mu')/ sum(likelihood .* mu');

    x_hat(:, 1) = imm(1).state.x;
    x_hat(:, 2) = imm(2).state.x;
    x_hat(:, 3) = imm(3).state.x;

    imm(1).state.k = k(1);
    imm(2).state.k = k(2);
    imm(3).state.k = k(3);

    x_out = x_hat * k';
    innov = x_hat - x_out;
    
    P_out = zeros(size(x_hat, 1));
    for j=1:n_mode
        P_out = P_out + k(j) * (imm(j).state.P + innov(:, j) * innov(:, j)');

        x_log{i}(:, j)    = imm(j).state.x;
        P_log{i}(:, :, j) = imm(j).state.P;
        k_log{i}(j)       = imm(j).state.k;

          
    end
   
    imm(1).x_out(:, i)    = x_out;
    imm(1).P_out(:, :, i) = P_out;
end



x_log_plot = zeros(4, n, n_mode);

for j = 1:n_mode
    for i=1:n
        x_log_plot(:, i, j) = x_log{i}(:, j);
    end
end

figure(); hold on; grid on;
colors = {'c', 'g', 'b'};
for j=1:n_mode
    plot(x_log_plot(1, :, j), x_log_plot(2, :, j), 'Color',colors{j}, "LineStyle","--"); hold on;
    legends{j} = (sprintf("imm filter %d", j));
end
legends{4} = "Ground Truth";
legends{5} = "IMM Mean";

plot(x(1, :), x(2, :), 'k', 'linewidth', 2);hold on;
plot(imm(1).x_out(1, :), imm(1).x_out(2, :), 'r*');
xlabel('x'); ylabel('y'); title('Trajectories');legend(legends);





