function [x] = genGauss(N, mu, scale, beta)
% % generate a genGauss vector of variates follows the pdf  f(x) = B exp( -|Ax|^beta )
% % The method is from https://blogs.sas.com/content/iml/2016/09/21/simulate-generalized-gaussian-sas.html
% % Inputs:
% %       N: dimension of the random vector
% %       mu: location of the GGN; real
% %       scale:  scale of the GGN; positive
% %       beta: shape of the GGN; positive

    % simulate for mu=0 and scale=1, then scale and translate
    a = sqrt(gamma(3/beta) / gamma(1/beta));
    b = a^beta;

    x = gamrnd(1/beta, 1/b, N, 1); % X ~ Gamma(1/beta, 1/b)
    sgn = binornd(1, 0.5, N, 1);
    sgn = -1 + 2*sgn; % random {-1, 1}

    x = mu + scale * sgn .* x.^(1/beta);
end



