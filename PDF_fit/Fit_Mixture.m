function S = Fit_Mixture(X)
% This code fits 1D epsilon data to a two-component lognormal mixture
% distribution in log10 space using maximum likelihood estimation (MLE)
% via the expectation-maximization (EM) algorithm.
%
% Input:
%   X            : 1D vector of dissipation data (epsilon), X > 0
%
% Output:
%   S            : struct containing data, empirical CDF, and fitted mixture
%                  distribution in log10-epsilon space
%
%   S.raw        : cleaned input data in linear epsilon space (after
%                  removing NaNs, non-positive values, and X >= 0.1)
%   S.ecdf       : empirical cumulative distribution function values of
%                  log10(S.raw)
%   S.x          : coordinates in log10-epsilon space
%                  corresponding to S.ecdf
%
%   S.w          : 1×K vector of mixture weights for each lognormal
%                  component in log10 space (sum(S.w) = 1)
%   S.mu         : 1×K vector of means of the Gaussian components in
%                  log10-epsilon space (sorted from smallest to largest)
%   S.sig        : 1×K vector of standard deviations of the Gaussian
%                  components in log10-epsilon space (sorted consistent
%                  with S.mu)
%
%   S.xfit       : row vector of grid points in log10-epsilon space on
%                  which the fitted mixture CDF and PDF are evaluated
%   S.cdf_mix    : fitted mixture CDF evaluated on S.xfit (monotonically
%                  non-decreasing, constrained to [0,1])
%   S.pdf_mix    : fitted mixture PDF in log10-epsilon space, evaluated on
%                  S.xfit(2:end) (obtained by numerical differentiation of
%                  S.cdf_mix)
%
%   S.pdf_comp   : K×(length(S.xfit)-1) array; each row k is the
%                  contribution of component k to the total PDF in
%                  log10-epsilon space, evaluated on S.xfit(2:end)
%
%   S.mean_mix   : scalar; fitted mean of epsilon in linear space
%                  (i.e. E[epsilon]) obtained by integrating S.pdf_mix over
%                  epsilon = 10^(log10 epsilon)



K = 2; %Two-component mixture, can be ajusted to higher ones if needed


%=======Step1: Modify Input======

%work in log10 space, convention for epsilons


X = X(~isnan(X) & X>0 & X<0.1); %Preclude suspicious data that are negative or too large
logX        = log10(X); 



%=======Step2: Lognormal Mixture Fit======

opt = statset('MaxIter',5e3,'TolFun',1e-10,'Display','final');
gm  = fitgmdist(logX,K,'RegularizationValue',1e-6,...
               'Replicates',10,'Options',opt);


% sort components by mean, so that first component has smallest mean
[mu_sorted, idx] = sort(gm.mu(:)');
w_sorted         = gm.ComponentProportion(idx);
sig_sorted       = sqrt(reshape(gm.Sigma,1,[])); 
sig_sorted = sig_sorted(idx);




%=======Step3: Prepare Output======


S.raw        = X;     %raw data;
[S.ecdf,S.x] = ecdf(logX); %Emperical cdf;


S.w   = w_sorted; %weights of two components
S.mu  = mu_sorted;
S.sig = sig_sorted;


dx     = 1e-5;
S.xfit = (min(S.x)-2):dx:(max(S.x)+2); %fitted pdf and cdf is defined on S.xfit grid
cdf_mix = zeros(size(S.xfit));
for k = 1:K
   cdf_mix = cdf_mix + S.w(k)*normcdf((S.xfit - S.mu(k))/S.sig(k));
end
cdf_mix       = min(max(cdf_mix,0),1);
S.cdf_mix     = cummax(cdf_mix);
S.pdf_mix     = diff(S.cdf_mix)/dx;



for k = 1:K
   S.pdf_comp(k,:) = S.w(k)*normpdf(S.xfit(2:end), S.mu(k), S.sig(k));
end
eps_grid     = 10.^(S.xfit(2:end));
S.mean_mix   = sum(S.pdf_mix .* eps_grid * dx); %Mean epsilon from the fitted distribution


end
