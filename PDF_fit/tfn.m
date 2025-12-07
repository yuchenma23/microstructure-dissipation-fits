function val = tfn(x, fx)
% TFN  Vectorized interface to Owen's T-function.
%
%   VAL = TFN(X, FX) computes Owen's T-function for each corresponding
%   element of X and FX, returning VAL of the same size. X and FX must be
%   the same size (or scalar‐expandable).
%
%   This wrapper calls an internal sub‐function TFN_SCALAR(...) for each
%   element. The sub‐function TFN_SCALAR is basically the original scalar
%   code from Young & Minder (Algorithm AS 76) / Burkardt, with slight
%   modifications for clarity.
%
% Licensing:
%   This code is distributed under the GNU LGPL license.
%
% Reference:
%   JC Young, Christoph Minder,
%   "Algorithm AS 76: An Algorithm Useful in Calculating Non-Central T and
%   Bivariate Normal Distributions," Applied Statistics, Vol. 23(3), 1974.

    % Ensure X and FX are compatible in size for elementwise operation.
    if ~isscalar(x) && ~isscalar(fx) && ~isequal(size(x), size(fx))
        error('TFN:SizeMismatch',...
              'Inputs X and FX must be the same size (or scalar).');
    end

    % Use arrayfun to apply the scalar TFN logic element‐by‐element.
    val = arrayfun(@tfn_scalar, x, fx);

end  % <-- End of the main (vectorized) TFN function


%-------------------------------------------------------------------------
function value = tfn_scalar(x, fx)
% TFN_SCALAR  Scalar version of Owen's T-function, for a single X, FX.
%
%    VALUE = TFN_SCALAR(X, FX) returns the scalar T-function of Owen.
%
%    Original references and licensing details are as in the main TFN
%    function. This sub‐function is meant to be called only via arrayfun
%    or similar vectorization wrappers.

  % Number of Gaussian quadrature points:
  ng = 5;

  % Gauss weights and abscissas:
  r = [ ...
      0.1477621, ...
      0.1346334, ...
      0.1095432, ...
      0.0747257, ...
      0.0333357 ];
  u = [ ...
      0.0744372, ...
      0.2166977, ...
      0.3397048, ...
      0.4325317, ...
      0.4869533 ];

  % Constants used in the algorithm:
  tp  = 0.159155;   % = 1/(2*pi)
  tv1 = 1.0E-35;
  tv2 = 15.0;
  tv3 = 15.0;
  tv4 = 1.0E-05;

  %-------------------------------------------------------------------------
  % 1) Test for X near zero:
  %-------------------------------------------------------------------------
  if abs(x) < tv1
    % value = (1/(2*pi)) * atan(fx)
    value = tp * atan(fx);
    return
  end

  %-------------------------------------------------------------------------
  % 2) Test for large values of |X|:
  %-------------------------------------------------------------------------
  if tv2 < abs(x)
    value = 0.0;
    return
  end

  %-------------------------------------------------------------------------
  % 3) Test for FX near zero:
  %-------------------------------------------------------------------------
  if abs(fx) < tv1
    value = 0.0;
    return
  end

  %-------------------------------------------------------------------------
  % 4) Possibly compute truncation point for FX if shape is large:
  %-------------------------------------------------------------------------
  xs = -0.5 * x * x;    % used in exponent
  x2 = fx;
  fxs = fx * fx;        % fx^2

  % If log(1 + fx^2) - (x^2/2)*fx^2 >= tv3, reduce x2 by Newton iteration:
  if tv3 <= ( log(1.0 + fxs) - xs*fxs )

    x1  = 0.5 * fx;
    fxs = 0.25 * fxs;

    while true

      rt = 1.0 + fxs;
      % The expression in parentheses is basically:
      %   (xs*fxs + tv3 - log(rt))
      % divided by (2*x1*(1/rt - xs))
      numerator   = xs*fxs + tv3 - log(rt);
      denominator = 2.0*x1*( (1.0/rt) - xs );

      x2_new = x1 + numerator / denominator;
      fxs_new = x2_new * x2_new;

      % Check convergence:
      if abs(x2_new - x1) < tv4
        x2  = x2_new;
        fxs = fxs_new;
        break
      end

      % Update for next iteration:
      x1  = x2_new;
      fxs = fxs_new;

    end
  end

  %-------------------------------------------------------------------------
  % 5) Gaussian quadrature on [0.5 - 0.5] with the formula from AS76:
  %-------------------------------------------------------------------------
  rt = 0.0;
  for i = 1:ng

    % Evaluate integrand at (0.5 + u(i)) and (0.5 - u(i)).
    arg1 = (0.5 + u(i));
    arg2 = (0.5 - u(i));

    r1 = 1.0 + fxs * arg1*arg1;
    r2 = 1.0 + fxs * arg2*arg2;

    rt = rt ...
        + r(i)* ( exp(xs * r1) / r1 + exp(xs * r2) / r2 );
  end

  value = rt * x2 * tp;

end  % <-- End of the scalar TFN sub‐function
