function piLambdaCond = lambdaFunc(lambda,i,beta,deltaS)
%LAMBDAFUNC Summary of this function goes here
%   Detailed explanation goes here
global alpha
global phi

piLambda = prod((beta.^(alpha)/gamma(alpha))*lambda.^(alpha-1)*exp(-beta.*lambda));
term1 = log(piLambda);
term2 = sum(sum(log(nbinpdf(deltaS,kappa,phi))));

piLambdaCond = term1 + term2;

end

