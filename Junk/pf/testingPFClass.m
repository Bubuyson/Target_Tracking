clc; clear; close all
% rng(5)

x0_M = [1000,1000,100,100];
P0 = diag([100,100,10,10])^2;

x = mvnrnd(x0_M,P0,1)';
x0 = x';
N=100;
T=1;
rng(6)

I2=eye(2);
Z2=zeros(2);
A = [I2 T*I2 ; Z2 I2];
B =[T^2/2*I2 ; T*I2];
Q=eye(2)*100;
R=diag([10^2,(25*pi/180)^2]); 

for i=2:N
   x(:,i) = A*x(:,i-1)+B*mvnrnd(zeros(1,2),Q,1)';
end

h=@(x) [sqrt(x(1,:).^2 + x(2,:).^2); atan2(x(2,:),x(1,:))]; 
y=h(x)+mvnrnd(zeros(1,2),R,N)';

f =@(T) [I2 T*I2 ; Z2 I2];
B =@(T) [T^2/2*I2 ; T*I2];

% x0 = [1000 1000 10 10];
P0 = diag([100 100 10 10].^2);
pf = ParticleFilter(2000, f, B, Q, R, h);
pf = pf.setResample("Multinomial", "MaxWeight", 5);
pf = pf.setApplyMH(true);
pf = pf.setApplyAuxiliary(true);
pf = pf.initialize(x0,P0);
tic
for i = 1:N
    pf = pf.update(y(:, i), 1);
end
toc
pf.plotTrajectory(x, [1, 2])
pf.calculateRMSE(x, [1 2; 3 4], 1)


