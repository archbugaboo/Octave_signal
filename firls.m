## Copyright (C) 2006 Quentin Spencer <qspencer@ieee.org>
##               2017 Ionescu Vlad <vlad.inf@gmail.com>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {@var{b} =} firls (@var{n}, @var{f}, @var{a})
## @deftypefnx {Function File} {@var{b} =} firls (@var{n}, @var{f}, @var{a}, @var{k})
## @deftypefnx {Function File} {@var{b} =} firls (@var{n}, @var{f}, @var{a}, @var{ftype})
## @deftypefnx {Function File} {@var{b} =} firls (@var{n}, @var{f}, @var{a}, @var{f2})
## @deftypefnx {Function File} {@var{b} =} firls (@var{n}, @var{f}, @var{a}, @var{k}, @var{ftype})
## @deftypefnx {Function File} {@var{b} =} firls (@var{n}, @var{f}, @var{a}, @var{ftype}, @var{f2})
##
## FIR filter design using least squares method. Returns a length @var{n}+1
## linear phase filter such that the integral of the weighted mean
## squared error in the specified bands is minimized.
##
## b = firls (@var{n}, @var{f}, @var{a}) creates a type I or II FIR with unity
## weights.
##
## b = firls (@var{n}, @var{f}, @var{a}, @var{k}) creates weighted types I and
## II FIRs.
##
## b = firls (@var{n}, @var{f}, @var{a}, @var{ftype}) creates all types of FIRs,
## differentiators, or Hilbert transformers, all with unity weights.
##
## b = firls (@var{n}, @var{f}, @var{a}, @var{f2}) creates types I and II FIRs
## with 1/f^2 weighting.
##
## b = firls (@var{n}, @var{f}, @var{a}, @var{k}, @var{ftype}) creates all
## types of FIRs, differentiators. Hilbert transformers can also be made, but
## the weights do not count.
##
## b = firls (@var{n}, @var{f}, @var{a}, @var{ftype}, @var{f2}) creates all
## types of FIRs, differentiators, and Hilbert transformers with 1/f^2
## weighting.
##
## The vector @var{f} specifies the frequencies of the band edges, normalized
## so that half the sample frequency is equal to 1.  Each band is specified by
## two frequencies, so the vector must have an even length.
##
## The vector @var{a} specifies the amplitude of the desired response at each
## band edge.
##
## The optional argument @var{k} is a weighting function that contains one
## value for each band that weights the mean squared error in that band. When
## not specified, it means the weights are one across all bands. It cannot be
## used together with @var{f2}.
##
## Also optional is @var{ftype}, of type string, which can be either null, '',
## or one of 'd', 'differentiator', 'h', 'Hilbert', or 'hilbert'. When
## specified, it will create a type III or IV FIR, which can also be a
## differentiator or a Hilbert transformer.
##
## The last optional argument is @var{f2}, which cannot be used together with
## @var{k}, specifies 1/f^2 weighting. It is the default option for
## differentiators, but not for the rest.
##
## @var{a} must be the same length as @var{f}, and @var{k} must be half the
## length of @var{f}. @var{n} can be odd or even, and the resulting filters do
## not increment the order behind the curtains in order to adjust the filter to
## a proper method. For example, a type II highpass is perfectly achievable,
## but the price to pay is the zero at Nyquist. The same goes for every other
## combination. It is up to the user to know what he or she wants, and how to
## achieve it.
##
## Examples with default (unity) weights:
##
## 30th order, type I highpass:
## h = firls(30, [0 0.3 0.4 1], [0 0 1 1]);
##
## 31st order, type II lowpass:
## h = firls(31, [0 0.3 0.4 1], [1 1 0 0]);
##
## 64th order, type III bandstop:
## h = firls(64, [0 0.3 0.4 0.7 0.8 1], [1 1 0 0 1 1]);
##
## 65st order, type IV bandpass:
## h = firls(51, [0 0.3 0.4 0.7 0.8 1], [0 0 1 1 0 0]);
##
## 48th order multiband:
## h = firls(48, [0 0.3 0.4 0.6 0.7 0.9], [0 1 0 0 0.5 0.5], 'h');
##
## Examples with weights:
##
## 42nd and 43rd order differentiators:
## h = firls(42, [0 0.3 0.4 1], [0 0.3 0 0]*pi, 'd');
## h = firls(43, [0 0.3 0.4 1], [0 0.3 0 0]*pi, [30 1], 'd');
##
## 49th and 50th order Hilbert transformers:
## h = firls(49, [0.1 1.0], [1 1], 'h'); # note the 1 at Nyquist
## h = firls(50, [0.1 0.9], [1 1], 10, 'h');
##
## Oddities:
##
## 26th order, type IV weighted lowpass:
## h = firls(26, [0 0.3 0.4 1], [1 1 0 0], [1 26], 'h'); # or 'd'
##
## 25th order, type II and IV highpass
## h = firls(25, [0 0.3 0.4 1], [0 0 1 1]);
## h = firls(25, [0 0.3 0.4 1], [0 0 1 1], 'h');
##
## 35th order, type II and type IV bandstop:
## h = firls(35, [0 0.3 0.4 0.7 0.8 1], [1 1 0 0 1 1], [1 10 1]);
## h = firls(35, [0 0.3 0.4 0.7 0.8 1], [1 1 0 0 1 1], [1 10 1], 'h');
##
## The least squares optimization algorithm for computing FIR filter
## coefficients is derived in detail in:
##
## I. Selesnick, "Linear-Phase FIR Filter Design by Least Squares,"
## http://cnx.org/content/m10577
## @end deftypefn

function h = firls(N, F, A, varargin);

# nr or arguments must be 3, 4, or 5
narginchk(3, 5);

# order must be a one one element vector...
if(length(N) != 1)
  error("The order (N) must be a vector of size 1.")
end
# ...and proper valued
if((N <= 0) || ischar(N))
  error("The order (N) must be positive definite.")
end

# handle the possible cases
oddN = mod(N, 2);

# three arguments => types I and II FIRs, unity weights
if(nargin == 3)
  K = ones(size(F(1:2:end)));
  fType = 0;
  f2 = 0;
end

# four arguments
if(length(varargin) == 1)
  if(ischar(varargin{1})) # diff or HT
    K = ones(size(F(1:2:end)));
    # check that ftype is a proper char
    switch(varargin{1})
      case {"h" "Hilbert" "hilbert"}
        fType = 1;
        f2 = 0;
      case {"d" "differentiator"}
        fType = 1;
        f2 = 1;
      case {"f2"}
        fType = 0;
        f2 = 1;
      otherwise
        print_usage()
    end
  else # type I or II FIR
    K = varargin{1};
    fType = 0;
    f2 = 0;
  end
end

# full spec
if(length(varargin) == 2)
  # make sure K is a number
  if(isnumeric(varargin{1}))
    K = varargin{1};
    if(strcmp(varargin{2}, 'h') || strcmp(varargin{2}, 'Hilbert') || ...
        strcmp(varargin{2}, 'hilbert') || strcmp(varargin{2}, 'd') || ...
        strcmp(varargin{2}, 'differentiator'))
      fType = 1;
    else
      error("'ftype' must be one of 'h', 'Hilbert', 'd', or 'differentiator'")
    end
    f2 = 0;
  else # check whether it's 'fType' or 'f2' (or not)
    if(strcmp(varargin{2}, 'h') || strcmp(varargin{2}, 'Hilbert') || ...
        strcmp(varargin{2}, 'hilbert') || strcmp(varargin{2}, 'd') || ...
        strcmp(varargin{2}, 'differentiator'))
      fType = 1;
      f2 = 1; # if it's 'fType', the 2nd can only be 'f2'
    else
      error("'ftype' must be one of 'h', 'Hilbert', 'd', or 'differentiator'")
    end
  end
end

# check the lengths of the vectors
if(length(F) != length(A))
  error("The sizes of the frequency and magnitude vectors must be equal.")
end

if((length(varargin) >= 1) && !ischar(varargin{1}) && ...
    (length(varargin{1}) != length(F)/2))
  error("The length of the weights vector must be half the length of the frequencies', or the amplitudes'.")
end

if(rem(length(F), 2) || rem(length(A), 2))
  error("The frequency and weight vectors must have an even size greater than or equal to 2.")
end

## Lengths alright? Check for correctness.

if(sum(diff(F) <= 0))
  error("The frequencies must be strictly increasing.")
end

if((F(1) < 0) || (F(end) > 1))
  error("The frequencies must lie in the interval [0..1], with 1 being Nyquist.")
end

if((length(F) != 2) && ((F(1) != 0) || (F(end) != 1)))
  error("The frequency vector must start at 0 and end with 1.")
end

# silently consider the integer part of N
N = fix(N);

# make sure the vectors are rows
A = A(:)';
F = F(:)';
## Make K the same length as F and A, to avoid indexing with floor((i+1)/2.
## Also make it alternating signs, to avoid the need of (-1)^n later on.
K = [-abs(K); abs(K)](:)';

# prepare a few helpers
M = floor(N/2);
bands = length(F);
w = F*pi;
A0 = fType*(1 - oddN);
pi2 = pi/2;
i1 = 1:2:bands;
i2 = 2:2:bands;
## Non-zero band check. Note: this doesn't work for Hilbert transformer.
bandTest = A(i1) + A(i2);
bandTest = [bandTest; bandTest](:)';

################################################################################
## Using 1/f^2 weighting means starting from the definition of q with W = 1/w^2
## and, with the help of Stegun's Handbook of Mathematical Functions, wxMaxima,
## and Wolframalpha:
##
##      ,-
##     /  cos(n*pi*f)                           cos(n*pi*f)
##    /  ------------- df = -n*pi*Si(n*pi*f) - ------------- + C
##   /        f*f                                    f
## -'
##
## where Si(x) is the sine integral. In Octave, expint() is defined, so:
##
##  Si(x) = imag(expint(1i*x)) + pi/2
##  Ci(x) = -real(expint(1i*x))
################################################################################
# calculate q
n = (0 : N+2*A0)'; # make it column vector from the start
q = zeros(size(n));
## Matlab uses 1/f^2 weighting only in the passband(s)
## TODO the two for loops can be merged at the cost of having for(if()), which
## whould be worse than if(for()), but which takes up more code, same for the
## b vector below; worth it?
ghostTweak = 4; # TODO empirically determined. <b>WHY?!</b>
## TODO if(for(if(if())))...?
if(f2) # 1/f^2 weighting
  for(m = 1:bands)
    if(bandTest(m)) # 1/f^2 weighting
      if(m == 1) # take care of F(1)=0
        q -= 0.0;
      else
        q -= ghostTweak*K(m)*(-pi*n.*(imag(E1(1i*n*w(m))) + pi2) - ...
          cos(n*w(m))/F(m));
      end
    else
      q -= K(m)*F(m)*sinc(n*F(m));
    end
  end
else
  for(m = 1:bands)
    q += K(m)*F(m)*sinc(n*F(m));
  end
end

# use q to build the Q matrix
Q = toeplitz(q(1 : M+1-A0)) + ...
  (-1)^fType*hankel(q(1+oddN+2*A0 : M+1+oddN+A0), q(M+1+oddN+A0 : N+1));

################################################################################
## In the same way q was derived, for the b vector, D(w) is a piecewise linear
## function, so it can be approximated as ax+b. The derivation is shown for
## sine, it's similar for cosine (needed only for types I and II):
##
##      ,-
##     /  (a*f + b)*sin(n*pi*f)
##    /  ----------------------- df = a*Si(n*pi*f) +
##   /              f*f
## -'                                    _                                 _
##                                      |                    sin(n*pi*f)   |
##                                    b*| pi*n*Ci(n*pi*f) - -------------  | + C
##                                      |_                       f        _|
##
## and ax+b is the linear interpolation of A and F vectors:
##
##  A[n+1] - A[n]
## ---------------*(x - F[n]) + A[n]
##  F[n+1] - F[n]
################################################################################
# adapted from original firls.m to include all types (I - IV) of filters.
n = (A0+0.5*oddN : M+0.5*oddN)'; # make column vector from the start
if(oddN)
  n2 = n;
  n3 = n;
else
  n2 = n(2:end);
  n3 = [1; n(2:end)];
end
## First check for 1/f^2, even if the order looks awkward, because the expint
## only needs to be calculated once, then picked apart with real() and imag().
## TODO if(for(if(if())))...?
if(f2) # 1/f^2 weighting
  sc = zeros(size(n));
  dif = zeros(size(n));
  tmp = zeros(size(n));
  for(m = 1:bands) # passband only
    l = m - 1 + mod(m, 2);
    if(bandTest(m)) # apply 1/f^2 weighting in the passband, only
      slope = (A(l+1) - A(l))/(F(l+1) - F(l));
      intercept = -slope*F(m) + A(m);
      if(m == 1) # take care of F(1)=0
        ## It seems there's no difference with/without tmp for F(1)=0, so just
        ## consider the singularities for Ci(x) and cos(x)/x as zero. This way
        ## numerical problems for large numbers are avoided in the case of
        ## e.g. Ci(1e-6), or cos(1e-6)/1e-6.
        sc += -K(m)*intercept*pi*n;
      else
        tmp = E1(1.0i*n*w(m));
        sc += K(m)*(slope*(imag(tmp) + pi2) + ...
          intercept*(pi*n.*(-real(tmp)) - sin(n*w(m))/F(m)));
      end
    else # stopband only
      AK = A(m)*K(m);
      if(!oddN && !fType) # odd lengths need special treatment
        sc += AK*[w(m); sin(n2*w(m))];
        dif += AK*[0.5*(w(l)^2 - w(l+1)^2)/(w(l+1) - w(l)); ...
          (cos(n2*w(l+1)) - cos(n2*w(l)))./(n2*(w(l+1) - w(l)))];
      else
        sc += AK*sin(n*w(m));
        dif += AK*(cos(n*w(l+1)) - cos(n*w(l)))./(n*(w(l+1) - w(l)));
      end
    end
  end
  # b vector
  b = ghostTweak*(dif + sc);
else # K decides the weighting
  sc = zeros(length(n), bands);
  dif = zeros(length(n), bands/2);
  if(fType) # types III, IV
    sc = -cos(n*w);
    dif = [sin(n*w(i2)) - sin(n*w(i1))]./(n*(w(i2) - w(i1)));
  else # types I, II
    if(!oddN) # odd lengths need special treatment
      sc = [w; sin(n2*w)];
      dif = [0.5*(w(i1).^2 - w(i2).^2)./(w(i2) - w(i1)); ...
        (cos(n2*w(i2)) - cos(n2*w(i1)))./(n2*(w(i2) - w(i1)))];
    else
      sc = sin(n*w);
      dif = (cos(n*w(i2)) - cos(n*w(i1)))./(n*(w(i2) - w(i1)));
    end
  end
  # b vector
  b = (kron(dif, [1, 1]) + sc)*(K.*A)(:)./(pi*n3);
end

# the rest of the algorithm
a = Q\b;

# form the impulse response
if(oddN)
  h = [(-1)^fType*flipud(a); a]';
else
  if fType
    h = [-flipud(a); 0; a]';
  else
    h = [a(end:-1:2); 2*a(1); a(2:end)]';
  end
end

endfunction

%% tests
%!error h = firls()
%!error h = firls(9)
%!error h = firls([1, 2])
%!error h = firls(9, 1)
%!error h = firls(9, 1, 2)
%!error h = firls(9, 1, 2, 3)
%!error h = firls(9, 1, 2, 3, 4)
%!error h = firls(9, 1, 2, 3, 4, 5)
%!error h = firls(9.9)
%!error h = firls(9, [])
%!error h = firls(9, [], [])
%!error h = firls(9, [], [], [])
%!error h = firls(9, [0 .2 .3 1], [1 2 3])
%!error h = firls(9, [.2 .5], 1)
%!error h = firls(9, 1, [1 2])
%!error h = firls(9, [-.1 .6 .9 1], [1 1 0 0])
%!error h = firls(-9, [0 .6 .9 1], [1 1 0 0])
%!error h = firls('x', [0 .6 .9 1], [1 1 0 0])
%!error h = firls(9, [0 .3 .6 1.3], [1 1 0 0])
%!error h = firls(9, [0 .6 .3 1], [1 1 0 0])
%!error h = firls(9, [0 .3 .6 1], [1 1 0 0], 1)
%!error h = firls(9, [0 .3 .6 1], [1 1 0 0], 'bla')
%!error h = firls("9", [0 .3 .6 1], [1 1 0 0])
