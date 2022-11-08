

%% Part 1
clear all, close all, clc
rng(5)
N = 6;

k = 3.5;
zn = [[0,0]',[k,0]',[0,k]',[0,-k]',[-k,0]'];
P = (1/20)*(diag(15*ones(1,5))+ones(5));

deltaT = 0.5;

alpha = 0.6;
sigmatil = [1 deltaT deltaT^2/2; 0 1 deltaT; 0 0 alpha];
sigmamat = [sigmatil, zeros(3);  zeros(3), sigmatil];

psiZtil = [deltaT^2/2 deltaT 0]';
psiWtil = [deltaT^2/2 deltaT 1]';
psiZ = [psiZtil zeros(3,1); zeros(3,1) psiZtil];
psiW = [psiWtil zeros(3,1); zeros(3,1) psiWtil];

X0 = mvnrnd(zeros(N,1), diag([500,5,5,200,5,5]))';

direction = round(rand(1)*4)+1;
Z0 = zn(:,direction);

sigma = 0.5;
Wi = mvnrnd(zeros(2,1), sigma^2*eye(2))';
Xi = X0;
Zi = Z0;

m = 3000;
states = zeros(m,2);

for i = 1:m
    Xi = sigmamat*Xi + psiZ*Zi + psiW*Wi;
    
    % update Zi and Wi
    update = rand(1);
    sum = 0;
    for j=1:5
        sum = P(direction,j) + sum;
        if update<sum
            direction = j;
            break
        end
    end
    Zi = zn(:,direction);
    Wi = mvnrnd(zeros(2,1), sigma^2*eye(2))';
    
    states(i,1) = Xi(1);
    states(i,2) = Xi(4);
    
end

plot(states(:,1), states(:,2));


%% Part 3
clear all, close all, clc

% ------- taken from problem 1 ---------

N = 100;
stat = 6;
k = 3.5;
zn = [[0,0]',[k,0]',[0,k]',[0,-k]',[-k,0]'];
P = (1/20)*(diag(15*ones(1,5))+ones(5));

deltaT = 0.5;

alpha = 0.6;
sigmatil = [1 deltaT deltaT^2/2; 0 1 deltaT; 0 0 alpha];
sigmamat = [sigmatil, zeros(3);  zeros(3), sigmatil];

psiZtil = [deltaT^2/2 deltaT 0]';
psiWtil = [deltaT^2/2 deltaT 1]';
psiZ = [psiZtil zeros(3,1); zeros(3,1) psiZtil];
psiW = [psiWtil zeros(3,1); zeros(3,1) psiWtil];

X0 = mvnrnd(zeros(stat,1), diag([500,5,5,200,5,5]),N)';

direction = round(rand(N,1)*4)+1;
Z0 = zn(:,direction);

sigma = 0.5;
Wi = mvnrnd(zeros(2,1), sigma^2*eye(2), N)';
Zi = Z0;

% ------- here starts SIS -------

load('RSSI-measurements.mat')
load('stations.mat')

v = 90;
eta = 3;
sigmaVn = 1.5;

% Vn = mvnrnd(0, sigmaVn);
% Y = v - 10*eta*log10(norm([X1, X2]' - pil)) + Vn;

m = size(Y,2)-1;

p = @(x1,x2,y) (1/sqrt(2*pi*sigmaVn^2))^N * prod(exp(-0.5*((y-v-10*eta*log10(norm([x1;x2]' - pos_vec)))/(sigmaVn^2))^2));
Xi = X0;
tau = zeros(2,m);
w = p(Xi(1,:), Xi(4,:), Y(:,1));

tau(:,1) = [Xi(1) Xi(4)].*w/w;
for k = 1:m-1
    Xi = sigmamat*Xi + psiZ*Zi + psiW*Wi;
    
    w = w*p(Xi(1),Xi(4),Y(k+1));
    tau(:,k+1) = [Xi(1);Xi(4)].*w./w;
    
    
    % update Zi and Wi
    update = rand(N,1);
    sum = 0;
    
    
    
    % Could move this outside the for-loop, do this before the time-loop
    for j=1:size(direction)
        
        sum = P(direction,j) + sum;
        if update<sum
            direction = j;
            break
        end
    end
    
    
    Zi = zn(:,direction);
    Wi = mvnrnd(zeros(2,1), sigma^2*eye(2),N)';

end





