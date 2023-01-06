# Jul 01 2022
I want to use the eSTA protocol to see if I can apply that with the josephson junction model.
The references are given in the paper 10.1103/PhysRevA.86.063623 where they have a full Hamiltonian and then they make some approximation in order to solve it.
The idea here is to evaluate the necessary corrections given by eSTA in order to obtain the solution to the full Hamiltonian.

- In the file `jj_calcs.nb` I will do those calculations to see if there is any semplification that arises
- I will then proceed evaluating the integrals that in turn will give me corrections.

I have found that the control parameter omega vanishes when I take the difference between the modified Hamiltonian and the STA one.
In theory I should be able to solve those integrals pretty easily.
Probably in the next days I will write down some latex notes.

# Jul 05 2022
After having discussed with me mates we decided to keep going down the route of eSTA so the idea is now to get the solutions of the STA protocol and evaluate the integrals that will give me the corrections to the modified Hamiltonian.
The approximated Hamiltonian they used in this paper is slighlty different from the actual harmonic oscillator but it is only a matter of constant factors.
I need to modify the eigenstates of the STA protocol accordingly and then solve the integrals.
I have tried to just substitute the mass with the  value h but I am not really happy with the result.
I could try to use some kind of change of variable to check if the STA eigentates have a better form.

# Jul 11 2022
I managed to solve the integrals and to obtain the Gns for `n` up to 10 and I found out that the only non zero terms are the ones of the form <2n|Hdiff|0>.
Now I need to evaluate the gradient of the fidelity by any means.
 
# Aug 02 2022
I still struggle to figure out if the assumption I made regarding the two levels and the left and right sites of the well is valid, so I will try to solve the harmonic oscillator Hamiltonian with the STA protocol and then I will try the very same protocol on the Hamiltonian where no approximations have been made.

# Aug 04 2022
I have been able to simulate the STA protocol for the quantum well, I also managed to write down the Hamiltonian in dimensionless units by changing the mass and the omega.
Now I need to transform the values they used in the paper to make sure the protocol works and with perfect fidelity.
In the paper they have ``\\omega_0 = 250 \\times 2\\pi Hz`` while `\\omega_f = 2.5 \\times 2\\pi Hz`. 
With these values of course ``\\gamma = 10``.
The total time of the protocol varies from 2ms up to 25 ms.
For the mass they used the one for Rb87.

# Aug 06 2022
I had the STA simulation running for the harmonic trap. I used the `Unitful.jl` package and made all the substitution.
I also defined the initial and final states by solving the time-independent Hamiltonian of the harmonic trap with the initial/final omega.
I used the ``timeevolution.schroedinger_dynamic`` function to time evolve the system and obtained almost perfect fidelity.
Next I need to set up the potential for the non approximated version of the Hamiltonian and time evolve that system with the STA protocol and check how much the fidelity dropped.

# Aug 08 2022
The solution for the harmonic oscillator cannot be used to solve the Hamiltonian for the bosonic system. I need to transform everything in terms of the number of particles ``N`` and ``\\Lambda``.
I need to do that first. It should not be super hard but I need to be careful.

I need to find a connection between the simulation of the working harmonic oscillator and the one for the Josephson junction. It's just a matter of numbers because the overall problem is the same.

# Aug 10 2022
I had the program running for both the two level picture and the angular momentum one, with agreeing results.
Now I need to express the STA protocol control paramters in such a way that I can carry on with the simulation.
I really need to write down some kind of notes where I work out some kind of change of variables because I really struggle with that.

# Aug 11 2022
Still struggling with the change of variable.
I will try to redo the calculation for the STA protocol only for the harmonic trap with single particle and try to mimic said calculation with the harmonic trap for the bosonic one.

# Aug 17 2022
I had the simulation running for both the Hamiltonians ``H_{BH}`` and ``H_{STA}``, the problem now is that I only have the eSTA correction for the Hamiltonian ``H_{N}`` obtained from the expansion of the initial Hamiltonian ``H_{BH}`` so I need to find out a way to either write down ``H_{BH}`` in terms of the position and momentum operators or try to solve ``H_{N}``.
Need to discuss that.

# Aug 22 2022	
I calculated the eSTA corrections `Gns` for a specific final time ``\neq 1 ``. Now I need to evaluate the other corrections `Kns` by solving the derivative of the non approximated Hamiltonian with respect to the control parameters.

I found out that there is a fast way to define the interpolating polynomial by using the Lagrange basis. I need to study that a little bit more but the results seem to be quite fast to evaluate.

In the Lagrange basis, if I want my polynomial to fit the points ``\lambda_{i}`` for ``i = 1,..., n``, I can write such polynomial as ``L(x) = \sum_{i=1}^{n} \lambda_i l_i(x)`` where ``l_i(x)`` are the so called *Lagrange basis* and are defined in a particular way. Those polynomials are produced in such a way that they are of the form ``\delta(x_i,x_j)`` i.e. they are not zero only for `` x = x_j`` so when I take the derivative of ``L(x)`` with respect of ``\lambda_j`` I am left with the corresponding ``l_j(x)``.

In order to calculate the gradient of ``H_{N}`` I have to add this polynomial correction ``L(t)`` to the control parameter ``\omega(t)`` and take the derivative with respect to each ``\lambda``. I will write the full calculations in the notes, but it is enough to say that the only surviving terms are the ones relative to ``z^2``.

# Oct 11 2022
In the following days I will try to find a way to connect the initial Bose Hubbard Hamiltonian and its approximated version.
My idea is to revert back the continuous variable z to  the discrete one that was used in the initial Hamiltonian, I will then need to multiply everything by ``JN``. I will probably find a solution to this problem if I check the paper where the approximations were carried out.

# Oct 17 2022
I managed to move from the Wikipedia formulation of the Bose Hubbard Hamiltonian to the one I got in the spin squeezed states paper.

# Nov 02 2022
I had some issues moving from the Bose Hubbard Hamiltonian to the STA one, so I will follow all the steps in the calculations to have some results that I am lately going to compare with the paper's ones. In particular I need to Taylor expand an exponential of a differential operator and apply that to the Taylor expansion of another function and I plan on doing that programmatically in mathematica. I will store all the progress I made in the `annotations.md` file.

# Nov 15 2022
Once I managed to obtain the continuous version of the Bose Hubbard Hamiltonian, I needed to evaluate Kns and Gns. 
The big problem in this case is that to calculate the Gns I could not obtain the analytical solution to the integral with respect of _z_ because of the `bh(z)` function. Luckily this was the only terms I had to integrate numerically, while the other two terms admit a closed form. 
Moreover for the Kns, I worked out a fast way to obtain the gradient so that I have to evaluate only one integral with respect of _z_ and it admit a closed form as well. 
I was lucky enough that I could also solve analytically the integrals without having to explicitely fix `m` so in theory I could have the code running even faster.
I have also tried to analytically solve the integrals with the `bh(z)` functions but with no luck. I need then to solve those integrals numerically.

Also, it looks like the two numeric integrals yield the same results. I need to find out if there is a problem in the code or if I can explain that. I did the difference between the integrals and they are not exactly zero, so the calculations are different but the results seem to be identical. I have to figure out why the two are the same.
# Nov 22 2022
Recently I moved to integrating with respect to the spatial coordinate over the whole real line instead of the interval [-1,1] as I found out the wavefunctions are normalized only if the domain of integration is the real line.
With this prescription, I found out that the integral <2|z^2|0> is the only one that is not zero and that simplifies a lot the calculations.

# Nov 29 2022
I tested the first round of corrections I obtained with the old version of the simulation program I wrote. The resulting fidelity is generally worse than the STA one, possibly because I used two different version of the constant factor _h_.
I will write down the calculations again to check where I can implement some improvements.

# Nov 30 2022
I redid the calculations and checked there was no mistakes and no matter what _h_ there should not be any problems as long as I am consistent with the definition.
I need to check how the integrals over the spatial variable change when I change the limits of integration. I could also use some formulas I found [ here ](https://giordano.github.io/Cuba.jl/stable/#Introduction) to change the variable of integration  between 0 and 1 instead of going all over the real line.

# Dec 08 2022
I started trying to evaluate the corrections with the new Hessian approach but unfortunately I obtain the corrections that are quite off. They are a hundred times bigger than the ones I obtain with the regular formula i.e. with no hessian involved.
It is probably due to some kind of constant factor that is missing in the calculations.

# Dec 13 2022   
The hessian seems to be off by a factor of 1/h^2. I do not know how to retrieve said results.

# Dec 29 2022
Recently I tried to use a simplified version of the gradient search for the eSTA protocol.
For each corrections I had, I multiplied each of them by a constant factor `\epsilon` and calculated the corresponding fidelity.
I then changed the value of `\epsilon` between 0.01 and 1.0 and checked how the fidelity landscape would change.
In general I found that the maximum fidelity was found for `\espilon` between 0.4 and 0.7, but the improvement was not relevant (max 10e-3).
I wrote the code and now I am in a position to do whatever I am asked.
