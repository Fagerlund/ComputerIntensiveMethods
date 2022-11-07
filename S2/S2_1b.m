
clear all, close all, clc


prob = @(x) exp((cos (x)^2 - 1));
trial = 1; 
accepted = false; 
while ~accepted
    Xcand = - pi/2 + pi*rand; 
    if rand < prob(Xcand)
        accepted = true; 
        X = Xcand;
    else 
        trial = trial + 1;
    end 
end 
