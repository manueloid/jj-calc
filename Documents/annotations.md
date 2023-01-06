# Jul 05 2022
I had some rough time trying to work out how to use the `Symbolics.jl` package. 
In particular I had some issues because some expressions of the form `sqrt(:symbol)^2` would not simplify correctly.
By searching online I ended up on an discourse page where they used the `Rewriter.jl` package and some of the features therein.
I have been trying to solve the problem fully in julia but the `Symbolics.jl` package as it seems like it is not fully implemented yet. I will probably move to mathematica or at least I will do some calculation by hand and then try to implement it.
I have tried to do all the calculation with full numerics but I am not happy with the result. I need to write in a nicer way the difference between H0 and Hn.

# Jul 07 2022
Since I immensely struggled to have `Symbolics.jl` working properly, I decided to move to `SymPy`.

Funnily enough there is a wrapper that works for julia as well so I am going to try to use that as it looks like SymPy is much more developed than `Symbolics.jl`. 

I found very interesting the non-evaluating derivatives `sympy.Derivatives` but very long to write and I found out you can actually alias a function in julia by writing something like

`const function = alias`

And it will work like a charm. For example in my case I just wrote

`const D = sympy.Derivative`
 
and I will save my future self a lot of time. I did not have to define some new function or anything, just used the *alias* capability of julia. 

# Jul 10 2022
I found out that it is usally better to use the `symbols` function when you need to define some variable that are identified by some kind of Greek letter as it makes it easier to print the latex form.
Moreover, you can use the format `@vars f()` to define a symbolic function f (instead of a variable)

# Jul 11 2022
I used the `.args` method of any expression variable in SymPy to split up the Hamiltonian as it was a summation of 5 terms.

Still, sympy works well but unfortunately not as well as mathematica when it comes to symbolics calculation.
I managed to have the mathematica file working alright but still it takes quite a while to analytically solve the integrals depending on `z` so I would assume it might be better to go full numerical from the get go (in such a way that I do not have to deal with mathematica and the plots will be easier to make)
I will think about it.

# Jul 14 2022
I have been playing around with `Symbolics` a little bit, here is what I found out:
1. If I use the imaginary Julia unit `im` and pair it with some kind of already defined variable inside one of those registered functions like `exp` or `sqrt` I will get an error.
2. If I use a dummy variable `I` for the imaginary unit and only then I make the substitution, I get the function to workand I also can work out the derivatives and stuff like that.

So here is what I will do from here on:
1. use the dummy variable `I` and then make the substitution `I=>im` straight after. 
2. Carry on with the calculation

This in theory should work, I will report back later in the day.
I had it working but I had to put in too much work and the results are not as fast as mathematica's

# Jul 15 2022	
I managed to have the integral calculation working and I obtained differnt gns.
I could try to speed it up a little bit but it is something I will work on if needs be.

# Aug 17 2022
I made some progress and now I am able to simulate the time evolution of the system for both the approximated version of the Josephon Junction Hamiltonian and the non approximated one.

I want to write down some code that will be able to return the fidelity of the STA protocol for different final times and then apply the very same control parameter to the non approximated Hamiltonian for the Josephon Junction and get the fidelity in said case.

I will split the discussion into two chunks to make it easier to read.

## Harmonic Oscillator
In this case the ingredients I mostly need are the solution for the Hamiltonian at initial and final time.
By changing the final time, the variables that will change are the control parameter and in turn the initial and final state.
The variables that will not change are the size of the box, the underlying basis and the representation of the momentum and position operators in different basis.
I will define a function that given the size of the box returns the two basis and then I will define the operators on said basis.
I would like to create a function that takes the boundary conditions and spits out the fidelity where the only changing parameter is the final time and I would like to achieve that without having to define the basis and the operator every time, since they do not change if I modify the final time.
 
To do that I think I have two options:
1. embed everything into a function where the initial part sets up the whole system, followed by a loop section where I evaluate the fidelity for tf in some kind of range
2. Create a function that takes only the final time as the input and produces the fidelity, assuming the operators have already been defined outside of the function scope. This is probably my choice because then I will be able to broadcast said function to a range of final time, speeding up the calculation. I need to remember to define the constant variables as `const`

## Josephon Junction Hamiltonian
I did that pretty easily and I saw the drop in fidelity 


## Note
I benchmarked the function control_ω and it was quite bad so I defined a new function that takes 1ns instead of 6ms and replaced it. Furthermore this new function control explicitely the final so it is no more set to 1, it changes every time.
Now I would like to benchmark the time evolution of the Schrodinger equation to see if I can obtain something faster.

For better performances it is always better to define the kind of input your function will take and also try to set up the variables as `const` whenever possible. Just for sake of the argument I have been able to reduce from 20ns to 1ns. 

# Aug 18 2022
Doing benchmark: of course the bigger the size of your array, the more memory you are going to use

I also changed the definition of `control_ω`

# Aug 19 2022
I found an error in the mathematica file so I found out that the calculation in fact takes a while to be carried out. 
I will hence try *again* to use julia.

# Aug 23 2022
No luck with the numeric evaluation, I will probably move back to Mathematica. However I found out I can change the time variable just by adding a `tf` factor in front of the derivative of `b(t)` in the definition of the STA wavefunction.

I switched to mathematica but it seems like I am still having some numerical bottlenecks, in particular when I consider a system with 100 particles. By reducing the number of particles, the instability vanishes.
I could try another change of variable but I do not feel like it.

# Aug 24 2022
Within the mathematica file I think I would have more luck if I solve the analytical integrals first and then assign the values, for both Gns and Zns.
There are too many numerical problems, I need to rewrite the Schrodinger equation of the bosonic oscillator and obtain the corresponding STA solutions.

# Aug 26 2022
I got numerical instabilities for N = 100 and ``\Lambda_f = 500`` when calculating the eSTA corrections. I probably will not be able to move everything into dimensionless units because there are some square roots in the Hamiltonian that are hardly solvable. I have then decided to scale everything back i.e. I set N =10 and ``Lambda_f = 50``. In this regime the STA protocol gets slightly worse as in the fidelity dropped from .9999999 to .995 but hopefully I will be able to get the corrections that will help me achieve better fidelity.

Don't know why but I got back to julia and decided to scale the time variable from 0 to 1.
1. The integral changes, and I need to add a factor tf in front of ``\int_0^t\omega_0/b(t)^2``
2. I need to include an extra 1/tf factor every time I take the derivative of b(t) with respect of time.

# Oct 06 2022
I found out how to include Latex code directly into Julia files. To do that  we need to append literal latex expression and `push!` them into the axis element

# Nov 22 2022
No big breakthrough until I started integrating over the whole real line for the spatial variable.
I then decided to start doing the calculations for the Kns and Gns in julia again.
It looks like it is working but I always have the feeling that my codes are too straightforward and can be optimized a lot more.
I need to find out how the hermite polynomials work 

# Nov 29 2022
I optimized as much as I could the `esta_calc.jl` file, trying  to use in place substitution and defining a new type holding all the relevant information of the protoctitution and defining a new type holding all the relevant information of the protocol.
Moreover I found the way to parallelize the calculation of some for loop via the `Threads.@threads` macro. It uses the number of cpu threads defined  by the environment variable `JULIA_NUM_THREADS` and automatically parallelizes the calculations. Very useful.
Now I need to check if I can parallelize the calculations of the integrals. That would mean the calculations will run much faster.

- Nice thing to know: creating N immutable struct variable is less costly than changing N times the value of one mutable type 

# Dec 01 2022
I have the program working quite fastly but the compilation time at the startup is quite huge. It looks like the compiler is spending a lot of time doing the so called _type inference_. 
This is the next thing I am going to work on.

# Dec 08 2022 
I am changing all the functions in such a way that they return values instead of anonymous functions, hopefully it will speed up the compilation time.
But now I have a problem while I try to broadcast an array onto a function that takes a custom type as  an input and it thus not work. I have to find out a solution for this feat.
-QuickFix- Use a for loop

