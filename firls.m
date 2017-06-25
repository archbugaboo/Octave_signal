## Copyright (C) 2006 Quentin Spencer <qspencer@ieee.org>
## Copyright (C) 2017 Ionescu Vlad <vlad.inf@gmail.com>
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
## @deftypefnx {Function File} {@var{b} =} firls (@var{n}, @var{f}, @var{a}, @var{k}, @var{ftype})
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
## b = firls (@var{n}, @var{f}, @var{a}, @var{k}, @var{ftype}) creates all
## types of FIRs, differentiators. Hilbert transformers can also be made, but
## the weights do not count.
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
## not specified, it means the weights are one across all bands.
##
## Also optional is @var{ftype}, of type string, which can be either null, '',
## or one of 'd', 'differentiator', 'h', 'Hilbert', or 'hilbert'. When
## specified, it will create a type III or IV FIR, which can also be a
## differentiator or a Hilbert transformer. The default option for
## differentiators is to have 1/f^2 weighting. If normal weighting is desired,
## use instead 'h'.
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
## h = firls (30, [0, 0.3, 0.4, 1], [0, 0, 1, 1]);
##
## 31st order, type II lowpass:
## h = firls (31, [0, 0.3, 0.4, 1], [1, 1, 0, 0]);
##
## 64th order, type III bandstop:
## h = firls (64, [0, 0.3, 0.4, 0.7, 0.8, 1], [1, 1, 0, 0, 1, 1]);
##
## 65st order, type IV bandpass:
## h = firls (65, [0, 0.3, 0.4, 0.7, 0.8, 1], [0, 0, 1, 1, 0, 0]);
##
## 48th order multiband:
## h = firls (48, [0, 0.3, 0.4, 0.6, 0.7, 0.9], [0, 1, 0, 0, 0.5, 0.5], 'h');
##
## Examples with weights:
##
## 42nd and 43rd order differentiators:
## h = firls (42, [0, 0.3, 0.4, 1], [0, 0.3, 0, 0]*pi, 'd');
## h = firls (43, [0, 0.3, 0.4, 1], [0, 0.3, 0, 0]*pi, [30, 1], 'd');
##
## 49th and 50th order Hilbert transformers:
## h = firls (49, [0.1, 1.0], [1, 1], 'h'); # note the 1 at Nyquist
## h = firls (50, [0.1, 0.9], [1, 1], 10, 'h');
##
## Oddities:
##
## 26th order, type IV weighted lowpass:
## h = firls (26, [0, 0.3, 0.4, 1], [1, 1, 0, 0], [1, 26], 'h'); # or 'd'
##
## 25th order, type II and IV highpass
## h = firls (25, [0, 0.3, 0.4, 1], [0, 0, 1, 1]);
## h = firls (25, [0, 0.3, 0.4, 1], [0, 0, 1, 1], 'h');
##
## 35th order, type II and type IV bandstop:
## h = firls (35, [0, 0.3, 0.4, 0.7, 0.8, 1], [1, 1, 0, 0, 1, 1], [1, 10, 1]);
## h = firls (35, [0, 0.3, 0.4, 0.7, 0.8, 1], [1, 1, 0, 0, 1, 1], [1, 10, 1], 'h');
##
## The least squares optimization algorithm for computing FIR filter
## coefficients is derived in detail in:
##
## I. Selesnick, "Linear-Phase FIR Filter Design by Least Squares,"
## http://cnx.org/content/m10577
## @end deftypefn

function h = firls (N, F, A, varargin);

## Nr or arguments must be 3, 4, or 5
narginchk (3, 5);

## Order must be a one one element vector...
if (length (N) != 1)
  error ("The order (N) must be a vector of size 1.")
endif
## ...and proper valued
if ((N <= 0) || ischar (N))
  error ("The order (N) must be positive definite.")
endif

## Handle the possible cases
oddN = mod (N, 2);

## Three arguments => types I and II FIRs, unity weights
if (nargin == 3)
  K = ones (size (F(1:2:end)));
  fType = 0;
  f2 = 0;
endif

## Four arguments
if (length (varargin) == 1)
  if (ischar (varargin{1})) # diff or HT
    K = ones (1, length (F)/2);
    fType = 1;
    switch (varargin{1})  # check that ftype is a proper char
      case {"h" "Hilbert" "hilbert"}
        f2 = 0;
      case {"d" "differentiator"}
        f2 = 1;
      otherwise
        print_usage()
    endswitch
  else  # type I or II FIR
    K = varargin{1};
    fType = 0;
    f2 = 0;
  endif
endif

## Full spec
if (length (varargin) == 2)
  fType = 1;  # it can only be one of the types III or IV
  if (isnumeric (varargin{1}))  # make sure K is a number
    K = varargin{1};
  else
    error ("The weights vector must be numeric.")
  endif
  switch (varargin{2})  # check that ftype is a proper char
    case {"h" "Hilbert" "hilbert"}
      f2 = 0;
    case {"d" "differentiator"}
      f2 = 1;
    otherwise
      print_usage()
  endswitch
endif

## Check the lengths of the vectors
if (length (F) != length (A))
  error ("The sizes of the frequency and magnitude vectors must be equal.")
endif

if ((length (varargin) >= 1) && ! ischar (varargin{1}) && ...
    (length (varargin{1}) != length (F)/2))
  error ("The length of the weights vector must be half the length of the frequencies', or the amplitudes'.")
endif

if (rem (length (F), 2) || rem (length (A), 2))
  error ("The frequency and weight vectors must have an even size greater than or equal to 2.")
endif

## Lengths alright? Check for correctness.
if (sum (diff (F) <= 0))
  error ("The frequencies must be strictly increasing.")
endif

if ((F(1) < 0) || (F(end) > 1))
  error ("The frequencies must lie in the interval [0..1], with 1 being Nyquist.")
endif

if ((length (F) != 2) && (F(1) != 0))
  error ("The frequency vector must start at 0.")
endif

if ((length (F) > 2) && mod (F, 2))
  error ("The frequency vector must have an even length.")
endif

N = fix (N);  # silently consider the integer part of N

## Make sure the vectors are columns
A = A(:)';
F = F(:)';
## Make K the same length as F and A, to avoid indexing with floor((i+1)/2.
## Also make it alternating signs, to avoid the need of (-1)^n later on.
K = [-abs(K); abs(K)](:)';

## Prepare a few helpers
M = floor (N/2);
bands = length (F);
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

## Calculate q
n = (0:N+2*A0)'; # make it column vector from the start
q = zeros (size (n));
## Matlab uses 1/f^2 weighting only in the passband(s)
## FIXME the two for loops can be merged at the cost of having for(if()), which
## whould be worse than if(for()), but which takes up more code, same for the
## b vector below; worth it?
ghostTweak = 4;  # FIXME empirically determined. <b>WHY?!</b>
## FIXME if (for (if (if ())))...?
if (f2)  # 1/f^2 weighting
  for (m = 1:bands)
    if (bandTest(m))  # 1/f^2 weighting
      if (m == 1)  # take care of F(1)=0
        q -= 0.0;
      else
        q -= ghostTweak*K(m)*(-pi*n.*(imag(E1(1i*n*w(m))) + pi2) - ...
          cos(n*w(m))/F(m));
      endif
    else
      q -= K(m)*F(m)*sinc(n*F(m));
    endif
  endfor
else
  for (m = 1:bands)
    q += K(m)*F(m)*sinc(n*F(m));
  endfor
endif

## Use q to build the Q matrix
Q = toeplitz(q(1:M+1-A0)) + ...
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

## Adapted from original firls.m to include all types (I - IV) of filters.
n = (A0+0.5*oddN : M+0.5*oddN)';  # make column vector from the start
if (oddN)
  n2 = n;
  n3 = n;
else
  n2 = n(2:end);
  n3 = [1; n(2:end)];
endif

## First check for 1/f^2, even if the order looks awkward, because the expint
## only needs to be calculated once, then picked apart with real() and imag().
## FIXME if (for (if (if ())))...?
if (f2)  # 1/f^2 weighting
  sc = zeros (size (n));
  dif = zeros (size (n));
  tmp = zeros (size (n));
  for (m = 1:bands)  # passband only
    l = m - 1 + mod (m, 2);
    if (bandTest(m))  # apply 1/f^2 weighting in the passband, only
      slope = (A(l+1) - A(l))/(F(l+1) - F(l));
      intercept = -slope*F(m) + A(m);
      if (m == 1)  # take care of F(1)=0
        ## It seems there's no difference with/without tmp for F(1)=0, so just
        ## consider the singularities for Ci(x) and cos(x)/x as zero. This way
        ## numerical problems for large numbers are avoided in the case of
        ## e.g. Ci(1e-6), or cos(1e-6)/1e-6.
        sc += -K(m)*intercept*pi*n;
      else
        tmp = E1(1.0i*n*w(m));
        sc += K(m)*(slope*(imag (tmp) + pi2) + ...
          intercept*(pi*n.*(-real (tmp)) - sin(n*w(m))/F(m)));
      endif
    else  # stopband only
      AK = A(m)*K(m);
      if (! oddN && ! fType)  # odd lengths need special treatment
        sc += AK*[w(m); sin(n2*w(m))];
        dif += AK*[0.5*(w(l)^2 - w(l+1)^2)/(w(l+1) - w(l)); ...
          (cos(n2*w(l+1)) - cos(n2*w(l)))./(n2*(w(l+1) - w(l)))];
      else
        sc += AK*sin(n*w(m));
        dif += AK*(cos(n*w(l+1)) - cos(n*w(l)))./(n*(w(l+1) - w(l)));
      endif
    endif
  endfor
  b = ghostTweak*(dif + sc);
else  # K decides the weighting
  sc = zeros (length (n), bands);
  dif = zeros (length (n), bands/2);
  if(fType)  # types III, IV
    sc = -cos(n*w);
    dif = [sin(n*w(i2)) - sin(n*w(i1))]./(n*(w(i2) - w(i1)));
  else  # types I, II
    if (! oddN)  # odd lengths need special treatment
      sc = [w; sin(n2*w)];
      dif = [0.5*(w(i1).^2 - w(i2).^2)./(w(i2) - w(i1)); ...
        (cos(n2*w(i2)) - cos(n2*w(i1)))./(n2*(w(i2) - w(i1)))];
    else
      sc = sin(n*w);
      dif = (cos(n*w(i2)) - cos(n*w(i1)))./(n*(w(i2) - w(i1)));
    endif
  endif
  b = (kron (dif, [1, 1]) + sc)*(K.*A)(:)./(pi*n3);
endif

## The rest of the algorithm
a = Q\b;

## Form the impulse response
if (oddN)
  h = [(-1)^fType*flipud(a); a]';
else
  if (fType)
    h = [-flipud(a); 0; a]';
  else
    h = [a(end:-1:2); 2*a(1); a(2:end)]';
  endif
endif

endfunction

%% expected output
%!test
%! x = [0.00469288328903675; ...
%!      -0.00163814821009947; ...
%!      -0.00981811774880307; ...
%!      -0.00827987877092834; ...
%!      0.00626749116017757; ...
%!      0.01958900144035464; ...
%!      0.01198122340365912; ...
%!      -0.01630409407768365; ...
%!      -0.03600919905667101; ...
%!      -0.01520878044033585; ...
%!      0.03919160528418460; ...
%!      0.07087566594952432; ...
%!      0.01740959934222501; ...
%!      -0.12559899212208905; ...
%!      -0.28304805248818982; ...
%!      0.64847607981652122; ...
%!      -0.28304805248818982; ...
%!      -0.12559899212208905; ...
%!      0.01740959934222501; ...
%!      0.07087566594952432; ...
%!      0.03919160528418460; ...
%!      -0.01520878044033585; ...
%!      -0.03600919905667101; ...
%!      -0.01630409407768365; ...
%!      0.01198122340365912; ...
%!      0.01958900144035464; ...
%!      0.00626749116017757; ...
%!      -0.00827987877092834; ...
%!      -0.00981811774880307; ...
%!      -0.00163814821009947; ...
%!      0.00469288328903675];
%! N = 30;
%! f = [0, 0.3, 0.4, 1];
%! A = [0, 0, 1, 1];
%! h = firls (N, f, A)';
%! assert (x, h, 1e-15);

%!test
%! x = [-5.72090574085339e-03; ...
%!      -2.46046647090348e-03; ...
%!      6.76837288733887e-03; ...
%!      1.14496771061725e-02; ...
%!      2.11533503057528e-03; ...
%!      -1.50529808933582e-02; ...
%!      -1.95870905692682e-02; ...
%!      9.88574373909403e-04; ...
%!      3.01841972156656e-02; ...
%!      3.18223810267019e-02; ...
%!      -1.05996614079784e-02; ...
%!      -6.27804784192314e-02; ...
%!      -5.71192981777838e-02; ...
%!      4.61676142678906e-02; ...
%!      2.09448476815995e-01; ...
%!      3.33494082295982e-01; ...
%!      3.33494082295982e-01; ...
%!      2.09448476815995e-01; ...
%!      4.61676142678906e-02; ...
%!      -5.71192981777838e-02; ...
%!      -6.27804784192314e-02; ...
%!      -1.05996614079784e-02; ...
%!      3.18223810267019e-02; ...
%!      3.01841972156656e-02; ...
%!      9.88574373909403e-04; ...
%!      -1.95870905692682e-02; ...
%!      -1.50529808933582e-02; ...
%!      2.11533503057528e-03; ...
%!      1.14496771061725e-02; ...
%!      6.76837288733887e-03; ...
%!      -2.46046647090348e-03; ...
%!      -5.72090574085339e-03];
%! N = 31;
%! f = [0, 0.3, 0.4, 1];
%! A = [1, 1, 0, 0];
%! h = firls (N, f, A)';
%! assert (x, h, 1e-15);

%!test
%! x = [-4.62543598284193e-04; ...
%!      7.08325591974923e-04; ...
%!      5.43663578714615e-05; ...
%!      7.17064531291712e-04; ...
%!      -1.32249184788049e-04; ...
%!      -2.92489836115204e-03; ...
%!      1.17008482578612e-03; ...
%!      5.55486543001852e-04; ...
%!      1.87446908977785e-03; ...
%!      3.69267617843608e-03; ...
%!      -7.34771373290465e-03; ...
%!      -1.85430762554092e-03; ...
%!      1.17415787975826e-03; ...
%!      2.80484430037048e-04; ...
%!      1.43705960494004e-02; ...
%!      -6.62280296853265e-03; ...
%!      -1.18708331189152e-02; ...
%!      1.12126501936955e-03; ...
%!      -9.90823255172805e-03; ...
%!      2.62512465483768e-02; ...
%!      1.34265772648076e-02; ...
%!      -2.66705673236500e-02; ...
%!      -8.29625890900374e-05; ...
%!      -3.24339022362673e-02; ...
%!      1.74694057721996e-02; ...
%!      7.06644242451211e-02; ...
%!      -3.34728491784523e-02; ...
%!      -1.78630854409270e-03; ...
%!      -7.02759026887394e-02; ...
%!      -9.19561342010886e-02; ...
%!      2.85134838675242e-01; ...
%!      6.05289017874539e-02; ...
%!      5.97403103359468e-01; ...
%!      6.05289017874539e-02; ...
%!      2.85134838675242e-01; ...
%!      -9.19561342010886e-02; ...
%!      -7.02759026887394e-02; ...
%!      -1.78630854409270e-03; ...
%!      -3.34728491784523e-02; ...
%!      7.06644242451211e-02; ...
%!      1.74694057721996e-02; ...
%!      -3.24339022362673e-02; ...
%!      -8.29625890900374e-05; ...
%!      -2.66705673236500e-02; ...
%!      1.34265772648076e-02; ...
%!      2.62512465483768e-02; ...
%!      -9.90823255172805e-03; ...
%!      1.12126501936955e-03; ...
%!      -1.18708331189152e-02; ...
%!      -6.62280296853265e-03; ...
%!      1.43705960494004e-02; ...
%!      2.80484430037048e-04; ...
%!      1.17415787975826e-03; ...
%!      -1.85430762554092e-03; ...
%!      -7.34771373290465e-03; ...
%!      3.69267617843608e-03; ...
%!      1.87446908977785e-03; ...
%!      5.55486543001852e-04; ...
%!      1.17008482578612e-03; ...
%!      -2.92489836115204e-03; ...
%!      -1.32249184788049e-04; ...
%!      7.17064531291712e-04; ...
%!      5.43663578714615e-05; ...
%!      7.08325591974923e-04; ...
%!      -4.62543598284193e-04];
%! N = 64;
%! f = [0, 0.3, 0.4, 0.7, 0.8, 1];
%! A = [1, 1, 0, 0, 1, 1];
%! h = firls (N, f, A)';
%! assert (x, h, 1e-15);

%!test
%! x = [8.49797987389339e-04; ...
%!      -4.78730552240640e-04; ...
%!      -5.26718382521342e-04; ...
%!      -2.58892977076172e-04; ...
%!      -1.20462250365382e-03; ...
%!      2.94299191816922e-03; ...
%!      1.01632146890130e-03; ...
%!      -1.90236309949536e-03; ...
%!      5.64526012032580e-05; ...
%!      -5.39327061898156e-03; ...
%!      2.80833818756474e-03; ...
%!      7.35759454705606e-03; ...
%!      -2.67440012053878e-03; ...
%!      2.44680106210254e-03; ...
%!      -1.01321091868340e-02; ...
%!      -7.57432004472130e-03; ...
%!      1.68183874260813e-02; ...
%!      8.32142385793243e-04; ...
%!      7.49916047426855e-03; ...
%!      -5.06890891793290e-03; ...
%!      -3.37108213479917e-02; ...
%!      1.72542434106304e-02; ...
%!      1.22888906832026e-02; ...
%!      1.17441016540351e-02; ...
%!      2.74252898348409e-02; ...
%!      -6.90324790184121e-02; ...
%!      -1.79621631153165e-02; ...
%!      3.47923045699895e-02; ...
%!      2.77925485139401e-03; ...
%!      1.38742303273865e-01; ...
%!      -9.50323503318502e-02; ...
%!      -2.92652299796673e-01; ...
%!      2.56069295992003e-01; ...
%!      2.56069295992003e-01; ...
%!      -2.92652299796673e-01; ...
%!      -9.50323503318502e-02; ...
%!      1.38742303273865e-01; ...
%!      2.77925485139401e-03; ...
%!      3.47923045699895e-02; ...
%!      -1.79621631153165e-02; ...
%!      -6.90324790184121e-02; ...
%!      2.74252898348409e-02; ...
%!      1.17441016540351e-02; ...
%!      1.22888906832026e-02; ...
%!      1.72542434106304e-02; ...
%!      -3.37108213479917e-02; ...
%!      -5.06890891793290e-03; ...
%!      7.49916047426855e-03; ...
%!      8.32142385793243e-04; ...
%!      1.68183874260813e-02; ...
%!      -7.57432004472130e-03; ...
%!      -1.01321091868340e-02; ...
%!      2.44680106210254e-03; ...
%!      -2.67440012053878e-03; ...
%!      7.35759454705606e-03; ...
%!      2.80833818756474e-03; ...
%!      -5.39327061898156e-03; ...
%!      5.64526012032580e-05; ...
%!      -1.90236309949536e-03; ...
%!      1.01632146890130e-03; ...
%!      2.94299191816922e-03; ...
%!      -1.20462250365382e-03; ...
%!      -2.58892977076172e-04; ...
%!      -5.26718382521342e-04; ...
%!      -4.78730552240640e-04; ...
%!      8.49797987389339e-04];
%! N = 65;
%! f = [0, 0.3, 0.4, 0.7, 0.8, 1];
%! A = [0, 0, 1, 1, 0, 0];
%! h = firls (N, f, A)';
%! assert (x, h, 1e-15);

%!test
%! x = [-0.001984287112619; ...
%!      0.005584703532062; ...
%!      -0.001099451383789; ...
%!      0.000296682062185; ...
%!      -0.006432628236318; ...
%!      0.001770872481965; ...
%!      -0.005537864709413; ...
%!      0.019419443183100; ...
%!      -0.003730039964359; ...
%!      -0.005016207687265; ...
%!      -0.016620436737913; ...
%!      0.004830573817572; ...
%!      -0.001356554505749; ...
%!      0.037294716882958; ...
%!      -0.000260971411954; ...
%!      -0.038031248400729; ...
%!      -0.021550384820893; ...
%!      0.002224553644658; ...
%!      0.036218442519860; ...
%!      0.060075475476243; ...
%!      0.034180066355077; ...
%!      -0.208651641643440; ...
%!      -0.058486399655920; ...
%!      -0.213732631805598; ...
%!      0.000000000000000; ...
%!      0.213732631805598; ...
%!      0.058486399655920; ...
%!      0.208651641643440; ...
%!      -0.034180066355077; ...
%!      -0.060075475476243; ...
%!      -0.036218442519860; ...
%!      -0.002224553644658; ...
%!      0.021550384820893; ...
%!      0.038031248400729; ...
%!      0.000260971411954; ...
%!      -0.037294716882958; ...
%!      0.001356554505749; ...
%!      -0.004830573817572; ...
%!      0.016620436737913; ...
%!      0.005016207687265; ...
%!      0.003730039964359; ...
%!      -0.019419443183100; ...
%!      0.005537864709413; ...
%!      -0.001770872481965; ...
%!      0.006432628236318; ...
%!      -0.000296682062185; ...
%!      0.001099451383789; ...
%!      -0.005584703532062; ...
%!      0.001984287112619];
%! N = 48;
%! f = [0, 0.3, 0.4, 0.6, 0.7, 0.9];
%! A = [0, 1, 0, 0, 0.5, 0.5];
%! h = firls (N, f, A, 'h')';
%! assert (x, h, 1e-15);

%!test
%! x = [-0.002501893069635; ...
%!      0.001621709433274; ...
%!      0.004608483193982; ...
%!      0.001079263379884; ...
%!      -0.006097531618565; ...
%!      -0.007282147564612; ...
%!      0.002024932544429; ...
%!      0.012524722603351; ...
%!      0.009404879378838; ...
%!      -0.008435528809980; ...
%!      -0.021678526274005; ...
%!      -0.009998524987725; ...
%!      0.020321329311149; ...
%!      0.035130738708294; ...
%!      0.007738646827657; ...
%!      -0.043218395455591; ...
%!      -0.060113093678052; ...
%!      -0.001715031879639; ...
%!      0.104553669943078; ...
%!      0.173028650212703; ...
%!      0.133392554052058; ...
%!      0.000000000000000; ...
%!      -0.133392554052058; ...
%!      -0.173028650212703; ...
%!      -0.104553669943078; ...
%!      0.001715031879639; ...
%!      0.060113093678052; ...
%!      0.043218395455591; ...
%!      -0.007738646827657; ...
%!      -0.035130738708294; ...
%!      -0.020321329311149; ...
%!      0.009998524987725; ...
%!      0.021678526274005; ...
%!      0.008435528809980; ...
%!      -0.009404879378838; ...
%!      -0.012524722603351; ...
%!      -0.002024932544429; ...
%!      0.007282147564612; ...
%!      0.006097531618565; ...
%!      -0.001079263379884; ...
%!      -0.004608483193982; ...
%!      -0.001621709433274; ...
%!      0.002501893069635];
%! N = 42;
%! f = [0, 0.3, 0.4, 1];
%! A = [0, 0.3, 0, 0]*pi;
%! h = firls (N, f, A, 'd')';
%! assert (x, h, 1e-15);

%!test
%! x = [2.14826965181594e-03; ...
%!      -5.78301863168259e-03; ...
%!      7.20640357680750e-04; ...
%!      7.07716530503224e-03; ...
%!      2.57492366522594e-03; ...
%!      -7.84620118500149e-03; ...
%!      -9.31356421086605e-03; ...
%!      3.32660462014073e-03; ...
%!      1.55374507608782e-02; ...
%!      8.84913283571647e-03; ...
%!      -1.34456072958950e-02; ...
%!      -2.42743909645316e-02; ...
%!      -3.83156892381569e-03; ...
%!      3.00740079188002e-02; ...
%!      3.42054375545578e-02; ...
%!      -9.43552581838081e-03; ...
%!      -6.01395600735065e-02; ...
%!      -4.98638670223835e-02; ...
%!      4.24990839563299e-02; ...
%!      1.52489257132067e-01; ...
%!      1.78690183372370e-01; ...
%!      7.94177383267942e-02; ...
%!      -7.94177383267942e-02; ...
%!      -1.78690183372370e-01; ...
%!      -1.52489257132067e-01; ...
%!      -4.24990839563299e-02; ...
%!      4.98638670223835e-02; ...
%!      6.01395600735065e-02; ...
%!      9.43552581838081e-03; ...
%!      -3.42054375545578e-02; ...
%!      -3.00740079188002e-02; ...
%!      3.83156892381569e-03; ...
%!      2.42743909645316e-02; ...
%!      1.34456072958950e-02; ...
%!      -8.84913283571647e-03; ...
%!      -1.55374507608782e-02; ...
%!      -3.32660462014073e-03; ...
%!      9.31356421086605e-03; ...
%!      7.84620118500149e-03; ...
%!      -2.57492366522594e-03; ...
%!      -7.07716530503224e-03; ...
%!      -7.20640357680750e-04; ...
%!      5.78301863168259e-03; ...
%!      -2.14826965181594e-03];
%! N = 43;
%! f = [0, 0.3, 0.4, 1];
%! A = [0, 0.3, 0, 0]*pi;
%! K = [30, 1];
%! h = firls (N, f, A, K, 'd')';
%! assert (x, h, 1e-15);

%!test
%! x = [-4.87357204240935e-05; ...
%!      -1.16485698858764e-04; ...
%!      -2.28339195699439e-04; ...
%!      -4.00746166178635e-04; ...
%!      -6.53519444635685e-04; ...
%!      -1.01006469923118e-03; ...
%!      -1.49761288252468e-03; ...
%!      -2.14749805900414e-03; ...
%!      -2.99555117011427e-03; ...
%!      -4.08272075476174e-03; ...
%!      -5.45609313793458e-03; ...
%!      -7.17058241788620e-03; ...
%!      -9.29172382131771e-03; ...
%!      -1.19002899931235e-02; ...
%!      -1.50999763728441e-02; ...
%!      -1.90304236621333e-02; ...
%!      -2.38899471263399e-02; ...
%!      -2.99769664161448e-02; ...
%!      -3.77701589956758e-02; ...
%!      -4.80964447361509e-02; ...
%!      -6.25231747456878e-02; ...
%!      -8.44225792300130e-02; ...
%!      -1.22590295286588e-01; ...
%!      -2.09336210716428e-01; ...
%!      -6.35657901226696e-01; ...
%!      6.35657901226696e-01; ...
%!      2.09336210716428e-01; ...
%!      1.22590295286588e-01; ...
%!      8.44225792300130e-02; ...
%!      6.25231747456878e-02; ...
%!      4.80964447361509e-02; ...
%!      3.77701589956758e-02; ...
%!      2.99769664161448e-02; ...
%!      2.38899471263399e-02; ...
%!      1.90304236621333e-02; ...
%!      1.50999763728441e-02; ...
%!      1.19002899931235e-02; ...
%!      9.29172382131771e-03; ...
%!      7.17058241788620e-03; ...
%!      5.45609313793458e-03; ...
%!      4.08272075476174e-03; ...
%!      2.99555117011427e-03; ...
%!      2.14749805900414e-03; ...
%!      1.49761288252468e-03; ...
%!      1.01006469923118e-03; ...
%!      6.53519444635685e-04; ...
%!      4.00746166178635e-04; ...
%!      2.28339195699439e-04; ...
%!      1.16485698858764e-04; ...
%!      4.87357204240935e-05];
%! N = 49;
%! f = [0.1, 1];
%! A = [1, 1];
%! h = firls (N, f, A, 'h')';
%! assert (x, h, 1e-15);

%!test
%! x = [-0.000133435071811; ...
%!      0.000000000000008; ...
%!      -0.000504709817584; ...
%!      0.000000000000020; ...
%!      -0.001338164675190; ...
%!      0.000000000000038; ...
%!      -0.002946651356234; ...
%!      0.000000000000059; ...
%!      -0.005752830944957; ...
%!      0.000000000000084; ...
%!      -0.010310633104538; ...
%!      0.000000000000107; ...
%!      -0.017355038244473; ...
%!      0.000000000000125; ...
%!      -0.027933173377369; ...
%!      0.000000000000134; ...
%!      -0.043759007605738; ...
%!      0.000000000000131; ...
%!      -0.068246589275604; ...
%!      0.000000000000114; ...
%!      -0.110110103028160; ...
%!      0.000000000000084; ...
%!      -0.201453570793047; ...
%!      0.000000000000045; ...
%!      -0.632962064594035; ...
%!      0.000000000000000; ...
%!      0.632962064594035; ...
%!      -0.000000000000045; ...
%!      0.201453570793047; ...
%!      -0.000000000000084; ...
%!      0.110110103028160; ...
%!      -0.000000000000114; ...
%!      0.068246589275604; ...
%!      -0.000000000000131; ...
%!      0.043759007605738; ...
%!      -0.000000000000134; ...
%!      0.027933173377369; ...
%!      -0.000000000000125; ...
%!      0.017355038244473; ...
%!      -0.000000000000107; ...
%!      0.010310633104538; ...
%!      -0.000000000000084; ...
%!      0.005752830944957; ...
%!      -0.000000000000059; ...
%!      0.002946651356234; ...
%!      -0.000000000000038; ...
%!      0.001338164675190; ...
%!      -0.000000000000020; ...
%!      0.000504709817584; ...
%!      -0.000000000000008; ...
%!      0.000133435071811];
%! N = 50;
%! f = [0.1, 0.9];
%! A = [1, 1];
%! K = 10;
%! h = firls (N, f, A, K, 'h')';
%! assert (x, h, 1e-15);

%!test
%! x = [-0.003022937893055; ...
%!      -0.010262389563569; ...
%!      -0.026905331240503; ...
%!      -0.048440190026946; ...
%!      -0.060457602098420; ...
%!      -0.049837464844988; ...
%!      -0.022306692637467; ...
%!      -0.007586962744112; ...
%!      -0.039285215295646; ...
%!      -0.120771780272181; ...
%!      -0.208208873470273; ...
%!      -0.233600887145252; ...
%!      -0.156107214739716; ...
%!      0.000000000000000; ...
%!      0.156107214739716; ...
%!      0.233600887145252; ...
%!      0.208208873470273; ...
%!      0.120771780272181; ...
%!      0.039285215295646; ...
%!      0.007586962744112; ...
%!      0.022306692637467; ...
%!      0.049837464844988; ...
%!      0.060457602098420; ...
%!      0.048440190026946; ...
%!      0.026905331240503; ...
%!      0.010262389563569; ...
%!      0.003022937893055];
%! N = 26;
%! f = [0, 0.3, 0.4, 1];
%! A = [1, 1, 0, 0];
%! K = [1, 26];
%! h = firls (N, f, A, K, 'h')';
%! assert (x, h, 1e-15);

%!test
%! x = [0.0178990170420163; ...
%!      -0.0309981382020000; ...
%!      0.0406264267627725; ...
%!      -0.0165946160117083; ...
%!      0.0387569103358096; ...
%!      -0.0684385398223223; ...
%!      0.0181960800188627; ...
%!      -0.0501928421261318; ...
%!      0.1306377719273110; ...
%!      -0.0331764282022474; ...
%!      0.0840065090785090; ...
%!      -0.4206296325450458; ...
%!      0.3010784434174334; ...
%!      0.3010784434174334; ...
%!      -0.4206296325450458; ...
%!      0.0840065090785090; ...
%!      -0.0331764282022474; ...
%!      0.1306377719273110; ...
%!      -0.0501928421261318; ...
%!      0.0181960800188627; ...
%!      -0.0684385398223223; ...
%!      0.0387569103358096; ...
%!      -0.0165946160117083; ...
%!      0.0406264267627725; ...
%!      -0.0309981382020000; ...
%!      0.0178990170420163];
%! N = 25;
%! f = [0 0.3 0.4 1];
%! A = [0, 0, 1, 1];
%! h = firls (N, f, A)';
%! assert (x, h, 1e-15);

%!test
%! x = [-0.00482858733565408; ...
%!      -0.01266780003649680; ...
%!      -0.00717045572499628; ...
%!      0.01187662485043596; ...
%!      0.02501268234043507; ...
%!      0.01041114060227501; ...
%!      -0.02681991742128989; ...
%!      -0.04762711603450899; ...
%!      -0.01298539054960517; ...
%!      0.06625082160141749; ...
%!      0.11317614947074235; ...
%!      0.01440992842973820; ...
%!      -0.54320408552794175; ...
%!      0.54320408552794175; ...
%!      -0.01440992842973820; ...
%!      -0.11317614947074235; ...
%!      -0.06625082160141749; ...
%!      0.01298539054960517; ...
%!      0.04762711603450899; ...
%!      0.02681991742128989; ...
%!      -0.01041114060227501; ...
%!      -0.02501268234043507; ...
%!      -0.01187662485043596; ...
%!      0.00717045572499628; ...
%!      0.01266780003649680; ...
%!      0.00482858733565408];
%! N = 25;
%! f = [0 0.3 0.4 1];
%! A = [0, 0, 1, 1];
%! h = firls (N, f, A, 'h')';
%! assert (x, h, 1e-15);

%!test
%! x = [-1.37292633719159e-02; ...
%!      6.54979171465630e-03; ...
%!      -2.69077740149984e-02; ...
%!      3.76306554359895e-02; ...
%!      -2.53489015186056e-02; ...
%!      4.56312464545050e-02; ...
%!      -4.00872361690424e-02; ...
%!      -1.89204036853308e-04; ...
%!      -2.33347427866400e-02; ...
%!      2.19335789807505e-02; ...
%!      1.27623656708041e-02; ...
%!      7.61230455487860e-02; ...
%!      -1.20184217197428e-01; ...
%!      7.10355471205113e-02; ...
%!      -2.06409230838411e-01; ...
%!      2.13451793662661e-01; ...
%!      9.97732350500253e-02; ...
%!      3.60851428052142e-01; ...
%!      3.60851428052142e-01; ...
%!      9.97732350500253e-02; ...
%!      2.13451793662661e-01; ...
%!      -2.06409230838411e-01; ...
%!      7.10355471205113e-02; ...
%!      -1.20184217197428e-01; ...
%!      7.61230455487860e-02; ...
%!      1.27623656708041e-02; ...
%!      2.19335789807505e-02; ...
%!      -2.33347427866400e-02; ...
%!      -1.89204036853308e-04; ...
%!      -4.00872361690424e-02; ...
%!      4.56312464545050e-02; ...
%!      -2.53489015186056e-02; ...
%!      3.76306554359895e-02; ...
%!      -2.69077740149984e-02; ...
%!      6.54979171465630e-03; ...
%!      -1.37292633719159e-02];
%! N = 35;
%! f = [0, 0.3, 0.4, 0.7, 0.8, 1];
%! A = [1, 1, 0, 0, 1, 1];
%! K = [1, 10, 1];
%! h = firls (N, f, A, K)';
%! assert (x, h, 1e-15);

%!test
%! x = [-1.25832412770046e-02; ...
%!      -1.77064540599929e-02; ...
%!      -3.92125618478823e-02; ...
%!      -2.26282755383741e-02; ...
%!      -2.44248934599686e-02; ...
%!      -1.16063929032390e-02; ...
%!      -3.21339192368622e-04; ...
%!      -5.39046658542452e-02; ...
%!      -4.25673740669957e-02; ...
%!      -6.31633052084354e-02; ...
%!      -5.93828386445542e-02; ...
%!      3.23769262586857e-02; ...
%!      -5.60377251242677e-02; ...
%!      -6.12674922029664e-02; ...
%!      -1.26303300495111e-01; ...
%!      -3.52124278520162e-01; ...
%!      -1.77941978680180e-02; ...
%!      -3.15612163936745e-01; ...
%!      3.15612163936745e-01; ...
%!      1.77941978680180e-02; ...
%!      3.52124278520162e-01; ...
%!      1.26303300495111e-01; ...
%!      6.12674922029664e-02; ...
%!      5.60377251242677e-02; ...
%!      -3.23769262586857e-02; ...
%!      5.93828386445542e-02; ...
%!      6.31633052084354e-02; ...
%!      4.25673740669957e-02; ...
%!      5.39046658542452e-02; ...
%!      3.21339192368622e-04; ...
%!      1.16063929032390e-02; ...
%!      2.44248934599686e-02; ...
%!      2.26282755383741e-02; ...
%!      3.92125618478823e-02; ...
%!      1.77064540599929e-02; ...
%!      1.25832412770046e-02];
%! N = 35;
%! f = [0, 0.3, 0.4, 0.7, 0.8, 1];
%! A = [1, 1, 0, 0, 1, 1];
%! K = [1, 10, 1];
%! h = firls (N, f, A, K, 'h')';
%! assert (x, h, 1e-15);

%% tests
%!error h = firls ()
%!error h = firls (9)
%!error h = firls ([1, 2])
%!error h = firls (9, 1)
%!error h = firls (9, 1, 2)
%!error h = firls (9, 1, 2, 3)
%!error h = firls (9, 1, 2, 3, 4)
%!error h = firls (9, 1, 2, 3, 4, 5)
%!error h = firls (9.9)
%!error h = firls (9, [])
%!error h = firls (9, [], [])
%!error h = firls (9, [], [], [])
%!error h = firls (9, [0 .2 .3 1], [1 2 3])
%!error h = firls (9, [.2 .5], 1)
%!error h = firls (9, 1, [1 2])
%!error h = firls (9, [-.1 .6 .9 1], [1 1 0 0])
%!error h = firls (-9, [0 .6 .9 1], [1 1 0 0])
%!error h = firls ('x', [0 .6 .9 1], [1 1 0 0])
%!error h = firls (9, [0 .3 .6 1.3], [1 1 0 0])
%!error h = firls (9, [0 .6 .3 1], [1 1 0 0])
%!error h = firls (9, [0 .3 .6 1], [1 1 0 0], 1)
%!error h = firls (9, [0 .3 .6 1], [1 1 0 0], 'bla')
%!error h = firls ("9", [0 .3 .6 1], [1 1 0 0],['a', 'b'])
%!error h = firls (9, [0 .6 .3 1], [1 1 0 0])