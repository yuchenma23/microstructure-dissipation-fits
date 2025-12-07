function F = lsn_cdf(z,p)
   alpha = p(1);  xi = p(2);  omega = p(3);
   zstar = (z - xi) / omega;
   % Skew-normal CDF:  Φ(z) − 2*T(z,α)
   F     = normcdf(zstar) - 2*tfn(zstar, alpha*ones(size(z)));
end
