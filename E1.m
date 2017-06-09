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
## @deftypefn  {Function File} {@var{y} =} expint (@var{x})
## @deftypefnx {Function File} {@var{y} =} expint (@var{x}, @var{eps})
##
## Exponential integral of the 1st order, according to Abramowitz&Stegun 5.1.1,
## defined as:
##
## @tex
## $$
## {\rm E_1} (x) = \int_x^\infty {e^{-t} \over t} dt
## $$
## @end tex
## @ifnottex
##
## @example
## @group
##            ,- inf
##           /
## E_1(x) = / exp(-t)/t dt
##         /
##       -' x
## @end group
## @end example
## @end ifnottex
##
## @var{x} can be real or complex. @var{eps} is the tolerance.
##
## The exponential integral is calculated with a recursion method, involving
## either a power series (with fast convergence for low values of @var{x}), and
## an continuous fraction expansion (with asymptotic convergence). The value of
## 4 has been chosen as a compromise, in terms of iterations, for swwitching
## between the two methods. If, after 100 iterations the result hasn't converged
## within the specified limits, the loop exits with the last obtained number.
##
## The code has been adapted from W. H. Press, S. A. Teukolsky, W. T. Vetterling
## "Numerical recipes. The Art of Scientific Computing".
##
## @end deftypefn

function y = expint(x, varargin)

narginchk(1, 2);

if((length(varargin) == 1) && ((ischar(varargin{1})) || (varargin{1} <= 0)))
  error('The tolerance must be a positive definite number.')
end

if(length(varargin) == 0)
  EPS = 1e-15;
else
  EPS = varargin{1};
end

y = zeros(size(x));
for(n = 1:numel(x))
  ## 4 seems a reasonable value for equally splitting the load
  if(abs(x(n)) > 4.0) # modfied Lentz splitting algorithm
    b = x(n) + 1.0;
    C = 1e30;
    D = 1.0/b;
    y(n) = D;
    for(k = 1:100)
      a = -(k*k);
      b += 2.0;
      D = 1.0/(b + a*D);
      C = b + a/C;
      Delta = C*D;
      y(n) *= Delta;
      if(abs(Delta - 1.0) <= EPS)
        #sprintf('%f,',x(n),n,k)
        y(n) *= exp(-x(n));
        break;
      end
    end
    #sprintf('The continuous fraction expansion failed after 100 iterations.')
    #y(n) *= exp(-x(n));
  else
    y(n) = -0.5772156649015328606065120900824-log(x(n));
    f = 1.0; # reusing same name variables
    for(k = 1:100)
      f *= -x(n)/k;
      C = -f/k;
      y(n) += C;
      if(abs(C) < abs(y)*EPS)
        #sprintf('%f,',x(n),n,k)
        #y(n) = y(n);
        break;
      end
    end
    #sprintf('The series expansion failed after 100 iterations.')
    #y(n) = y(n);
  end
end

endfunction
