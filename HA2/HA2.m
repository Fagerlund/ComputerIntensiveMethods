%% Constants and instantiating
clear all; close all; clc

iran_inf = readtable("iran_infected.csv"); iran_rec=readtable("iran_removed.csv");
germany_inf = readtable("germany_infected.csv"); germany_rec=readtable("germany_removed.csv");
germany_population = 83e6; iran_population = 84e6;
iran_inf = iran_inf{:,:}; iran_rec = iran_rec{:,:};
germany_inf = germany_inf{:,:}; germany_rec = germany_rec{:,:};

iran_sus = ones(size(iran_inf,1),1)*iran_population;
iran_sus = iran_sus - iran_inf - iran_rec;

germany_sus = ones(size(germany_inf,1),1)*germany_population;
germany_sus = germany_sus - germany_inf - germany_rec;

% Set initial and last value
iran_rec = [0; iran_rec]; 
iran_sus = [iran_population-iran_inf(1); iran_sus];
iran_inf = [iran_inf; 0];

germany_rec = [0; germany_rec]; 
germany_sus = [germany_population-germany_inf(1); germany_sus];
germany_inf = [germany_inf; 0];

% Constants
global phi
global alpha
phi = 0.995;
alpha = 2;


%% Fucntions



%% Iran

deltaI = diff(iran_inf);
deltaR = diff(iran_rec);
deltaS = -diff(iran_sus);
sum1 = 0; sum2 = 0;

% Parameters for investigation
breaks = 3;
beta = 2*ones(breaks,1);
a = 2; b = 3;
sigma = 5; M = 3;

% Initiate lambda, t, and pir
lambda = ones(breaks,1);
t = linspace(0,T,breaks+1)';
t = t(2:end-1);
pir = zeros(length(iran_inf),1);

sumx = 0;
sumit = 0;

iter = 1000;

for k = 1:iter
    % Metrolpolis-hastings for lambda
    % Create candidates for all lambda
    lambdaCand = lambda(:,k) + sigma*normrnd(0,1,breaks,1);
    
    for i = 1:breaks 
        alphaLambda = min(1,lambdaFunc(lambdaCand,i,beta,deltaS)./lambdaFunc(lambda(:,k),i,beta,deltaS));
        % if-statement
        if unifrnd(0,1) < alphaLambda(i)
            lambda(i,k+1) = lambdaCand(i);
        else
            lambda(i,k+1) = lambda(i,k);
        end
    end
    
    
    % Metrolpolis-hastings for t
    % Create candidates for all t
    tCand = t(:,k) + randi([-M,M],breaks-1,1);
    
    for i = 1:breaks-1
        alphaT = min(1,tFunc(tCand,i)./tFunc(t(:,k),i));
        % if-statement
        if unifrnd(0,1) < alphaT(i)
            t(i,k+1) = tCand(i);
        else
            t(i,k+1) = t(i,k);
        end
    end
    
    % Calculate pir
    betaA = sum(deltaS+deltaI)+a;
    betaB = sum(iran_inf(0:end-1)-deltaS-deltaI)+b;
    pir(k+1) = betarnd(betaA,betaB);
    
end


% for k = 1:height(Iran_inf)-1
%     epsilon1=normrnd(0,1); epsilon2=randi([-3,3]) ;
%     x=deltaI+deltaS;
%     sumx=sumx+x;
%     sumit=sumit+Iran_inf(k,1);
%     pir=betarnd(sumx-a,sumit-sumx+b); %Gibbs sampler,?? Or do we sample x too?
%     pir(k,1)=pir;
%     
%     %Random walk
%     %alpha=min(1,f(xstar)/f(xk));
%     
% end