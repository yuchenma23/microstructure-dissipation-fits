function S = Fit_LSN(X)
%This code fit a 1D epsilon data into the log-skew-normal distribution
%using the most-likelihood estimate method (MLE)
%
%Input:  X is a 1D vector of dissipation data (epsilon), X > 0
%
%Output: S is a struct containing data, empirical CDF, and fitted LSN curves
%   S.raw      : original input data X (in linear epsilon space)
%   S.ecdf     : empirical cumulative distribution function values of log10(X)
%   S.x        : support points (coordinates) in log10-epsilon space for S.ecdf
%   S.xfit     : dense grid in log10-epsilon space used to evaluate the fitted LSN CDF/PDF
%   S.cdf_fit  : fitted log-skew-normal CDF evaluated on S.xfit (log10-epsilon space)
%   S.pdf_fit  : fitted log-skew-normal PDF evaluated on S.xfit(2:end) (log10-epsilon space)
%   S.mean_fit : fitted mean of epsilon in linear space, i.e. E[epsilon] from the LSN fit


%=======Step 1: Modify Input======


X = X(~isnan(X) & X>0 & X<0.1); %Preclude suspicious data that are negative or too large
logX   = log10(X); %work in log10 space, convention for epsilons


%========Step 2: LSN Fit======
% initial guess of parameters from data
mu0    = mean(logX); 
s0     = std(logX);
if ~isfinite(s0) || s0<=0, s0 = 1; end
sk0    = mean((logX - mu0).^3)/(s0^3 + eps);  % sample skewness
a0     = max(-8, min(8, 2*sk0));            % mild map to alpha
q0     = [a0, mu0, log(max(s0, 1e-3))];     % [alpha0, xi0, omega0]


% use Nelder-Mead simplex algorithm to search for the best fit parameter
tiny = 1e-300;

negloglik = @(q) -sum( log( max(tiny, (2./exp(q(3))) .* ...
                 normpdf((logX - q(2))./exp(q(3))) .* ...
                 normcdf(q(1) .* (logX - q(2))./exp(q(3))) )) ); 
% pdf(x | alpha, xi, omega) = (2/omega) * phi((z-xi)/omega) * Phi(alpha*(z-xi)/omega)

opts = optimset('Display','iter','TolX',1e-8,'TolFun',1e-8, ...
               'MaxIter',5e4,'MaxFunEvals',5e4);
best_fit = fminsearch(negloglik, q0, opts);

alpha = best_fit(1);
xi    = best_fit(2);
omega = exp(best_fit(3));               



%======Step 3: Prepare Output======

S.raw      = X;
[S.ecdf, S.x] = ecdf(logX);   %coordinate for ecdf


dx         = 1e-5;
S.xfit     = (min(S.x)-2):dx:(max(S.x)+2); %coordinate for pdf/cdf of fitted curve

S.cdf_fit  = lsn_cdf(S.xfit,[alpha xi omega]);  % CDF in z-space
S.cdf_fit  = cummax( S.cdf_fit );               % force non-decreasing
S.pdf_fit  = diff(S.cdf_fit)/dx;                % PDF in z-space, simplified way of calculation, error is tiny


% mean in linear space
eps_fit    = 10.^(S.xfit(2:end));
S.mean_fit = sum(S.pdf_fit .* eps_fit * dx);


fprintf('[MLE] alpha = %.4f, xi = %.4f, omega = %.4f,  mean_fit = %.4e\n', ...
       alpha, xi, omega, S.mean_fit);

end