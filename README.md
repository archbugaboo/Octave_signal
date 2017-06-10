# Octave_signal
Some changes to Octave's signal package

Should have made this before, but now, here's the first github'ed version of Octave/signal/firls.m.
The name is lsfir2.m to be able to be used separately of firls.m, e.g. lsfir2(11,[0 .3 .4 1],[1 1 0 0],[1 10]).
When/if the script is accepted, it can be renamed afterwards.

There numbers should be spot on, compared to the Matlab script, for everything except the 1/f^2 weighting (differentiator). Two such examples have been provided in the mailing list:

>>> h = firls(11, [0 200 300 512]/512, [0 1 0 0], 'd');h'
> ans =
>    0.025308869942737
>   -0.018174108886369
>   -0.083194157264830
>    0.063294222003005
>    0.270747839846565
>    0.167064335397056
>   -0.167064335397056
>   -0.270747839846565
>   -0.063294222003005
>    0.083194157264830
>    0.018174108886369
>   -0.025308869942737

and this is my result:
h =

   0.0255707581685783
  -0.0183683660855015
  -0.0831834477531340
   0.0634034245580697
   0.2707190430504273
   0.1670027964290154
  -0.1670027964290154
  -0.2707190430504273
  -0.0634034245580697
   0.0831834477531340
   0.0183683660855015
  -0.0255707581685783
  
and:

>>> h = firls(10, [0 200 300 512]/512, [0 1 0 0], 'd');h'
> ans =
>    0.035928864351330
>   -0.067308010314644
>   -0.050456898575185
>    0.185894253348846
>    0.280660645384262
>                    0
>   -0.280660645384262
>   -0.185894253348846
>    0.050456898575185
>    0.067308010314644
>   -0.035928864351330

h =

   0.036180025458648
  -0.067401855849468
  -0.050570602450220
   0.185935657337846
   0.280795635112621
   0.000000000000000
  -0.280795635112621
  -0.185935657337846
   0.050570602450220
   0.067401855849468
  -0.036180025458648
  
This happens, most probably, because of ways to circumvent Ci(0) (cosine integral) and cos(0)/0. How, I don't know, but, initially, I forced the starting frequency to be 1e-6 instead of zero, but that yielded large numbers for the q and b vectors, which, no doubt, took their toll on precision. Now I simply forced the values to come out as zero for cos(0)/0, and Ci(0) for the cosine integral.

---

Octave's expint() has numerical problems for purely imaginary arguments (which are needed here), and they go worse as the argument's value increases. Temporary solution: create a separate variant for the exponential integral, E1(x), and use that.

With this change, the numbers come out just as with wxMaxima's and, given the existent examples, I dare say my results are better than Matlab's. If I am wrong, I am wrong, my apologies to the giant, but I only say this after I plot diff(abs(fft(h))), which shows a clear difference between the two results, with a much better ripple towards DC, and a visible 1/f^2 progression of the ripples.
