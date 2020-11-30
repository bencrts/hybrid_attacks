
# Hybrid Attacks on Small-secret LWE

This repository contains code to estimate the complexity of the Hybrid Decoding [How07] and Hybrid Dual [Alb17,CHHS19] attacks on small-secret LWE. This code is adapted from source code developed for the publications [APS15,ACW19,CP19], to which all authors have contributed -- however, all subsequent mistakes are mine. Pull requests are welcome -- if you find a bug/error feel free to report an issue. Note that for large parameters the Hybrid Decoding estimates can be fairly memory expensive, we are working on improving both the time and memory usage of this code for future releases. After cloning this repository, be sure to update all submodules:

`` git submodule update --init --recursive ``

to ensure that the code works correctly.

In this README we give a short overview of the Hybrid Decoding and Hybrid Dual attacks, and refer to the relevant literature for further details. We have orgainsed this README into separate sections for the Hybrid Decoding and Hybrid Dual attacks. In each section we explain the assumptions under which the estimates are made and how the code retrieves estimates. We also give an example usage of the code.

[Hybrid Decoding](#the-hybrid-decoding-attack) <br>
[Hybrid Dual](#the-hybrid-dual-attack)

There are several other implementations of hybrid attack estimates and we provide links to these alongside the bibliography. 
[Alternative Code and References](#alternative-code-and-references) 

# The Hybrid Decoding Attack

Given <img src = "images/m.svg" class = "inline"> LWE samples of dimension <img src = "images/n.svg" class = "inline">, we begin by constructing the lattice with basis matrix 

<img src = "images/basis_matrix.svg" class = "inline">

Here <img src = "images/lwe_matrix_cols_dropped.svg" class = "inline"> represents the public LWE matrix 
<img src = "images/lwe_matrix.svg" class = "inline"> with the first <img src = "images/tau.svg" class = "inline"> columns removed. If the first <img src = "images/tau.svg" class = "inline"> components of the LWE secret are zero, then the lattice <img src = "images/L(B).svg" class = "inline"> contains the vector <img src = "images/s_As.svg" class = "inline">, which is seperated from <img src = "images/zero_b.svg" class = "inline"> by the unusually short vector <img src = "images/s_minus_e.svg" class = "inline">. We can therefore find the vector <img src = "images/s_minus_e.svg" class = "inline"> by solving the Bounded Distance Decoding problem (BDD) around <img src = "images/zero_b.svg" class = "inline">. 

Otherwise, we have to begin a guessing procedure over the first <img src = "images/tau.svg" class = "inline"> components of the secret, which, for example, in the case of a uniformly random ternary secret would be <img src = "images/ternary_search_space.svg" class = "inline">. For each guess <img src = "images/vg_in_S.svg" class = "inline">, where <img src = "images/guess.svg" class = "inline">, and <img src = "images/unit_vector.svg" class = "inline"> represents the kth unit vector and we have <img src = "images/coefficients.svg" class = "inline">, then we solve BDD on the input point <img src = "images/BDD_input.svg" class = "inline"> in an attempt to find <img src = "images/BDD_output.svg" class = "inline">.

Thus, the Hybrid Decoding attack contains three phases:

1. Fixing the guessing dimension <img src = "images/tau.svg" class = "inline"> as above.
2. Performing lattice reduction (i.e. BKZ) on the lattice <img src = "images/L(B).svg" class = "inline">.
3. Solving a batch of BDD instances, using Babai's Nearest Plane algorithm, on input points derived from the target guess as described above.

This guessing phase is typically realised using a meet-in-the-middle process, which requires exponential memory, but allows for the searching phase to be performed in a more (time) efficient manner. To simulate this meet-in-the-middle process, we conservatively assume a square-root speed-up, and ignore any probability loss introduced. Thus, the estimates provided by our script are explicit underestimates of security, based on the assumptions outlined below. 

## Assumptions

The assumptions in our code are conservative. A discussion on assumptions used in the Hybrid Decoding attack can be found in [Wun16, ACW19, CP19].

| Assumption | Choice | Notes |
|---|---|---|
| Reduced basis shape | Geometric Series Assumption |   |
| Cost of Babai's algorithm |  <img src = "images/babai_cost.svg" class = "inline"> |   |
| Success probability of Babai's algorithm |  <img src = "images/babai_prob.svg" class = "inline"> |  See [Wun19]. Here <img src = "images/euler_beta.svg" class = "inline"> is the Beta function and <img src = "images/r_i.svg" class = "inline">.  |
| Meet in the middle speed-up | square-root |   |
| Meet in the middle success probability | 1 |  |
| Whether to consider cost of memory | No |  |

## What the Code Does

There are two important functions inside hybrid_decoding.py. The first is a function `hybrid_decoding_attack()` which can be used to generate the cost of the Hybrid Decoding attack on a given parameter set for a specified value of <img src = "images/tau.svg" class = "inline"> and <img src = "images/beta.svg" class = "inline">. 

The second is the function `parameter_search()`, which optimizes the complexity of the Hybrid Decoding attack over the choice of parameters <img src = "images/tau.svg" class = "inline"> and <img src = "images/beta.svg" class = "inline"> in the following way. We search over `range(60, n, 50)` for a value `beta'` and over `range(0, n, n//10)` for a value `tau'`. After this "optimal" pair `(tau', beta')` has been found, we search again over `range(beta' - 25, beta' + 25, 10)` for an optimal `beta''` and `range(tau' - n//20, tau'+n//20, n//100)` for an optimal `tau''`. The final output estimate is the function `hybrid_decoding_attack()` called with `tau = tau''` and `beta = beta''`.  There is a definitely a balance to be struck here as we want to search the parameter space reasonably well. However, for values of n seen in homomorphic encryption schemes (sometimes >= 2^16), it can be very costly to perform a granular search (particularly for the parameter <img src = "images/tau.svg" class = "inline">). Instructions to change the way this search is performed can be found in the source code, if a more granular search is desired.

## Example usage of the code

To retrieve a complexity estimate for the hybrid attack with blocksize <img src = "images/beta.svg" class = "inline"> = 100 and guessing dimension <img src = "images/tau.svg" class = "inline"> = 250 assuming we estimate lattice reduction using the BKZ.sieve cost model from [Est20]. This code can be run in the SageMathCell here.

```
%attach example_params.py
%attach hybrid_decoding.py
hybrid_decoding_attack(beta = 100, tau = 250, mitm = True, reduction_cost_model = est.BKZ.sieve, **example_64())

         rop:   2^65.1
         pre:   2^64.8
        enum:   2^62.5
        beta:      100
         |S|:   2^73.1
        prob: 0.104533
       scale:   12.760
          pp:       11
           d:     1798
      repeat:       42
```

To retrieve an estimate with optimised attack parameters (over the choice of <img src = "images/tau.svg" class = "inline"> and <img src = "images/beta.svg" class = "inline">) we call the function parameter_search().

```
%attach example_params.py
%attach hybrid_decoding.py
parameter_search(mitm = True, reduction_cost_model = est.BKZ.sieve, **example_64())

         rop:   2^63.4
         pre:   2^63.0
        enum:   2^61.5
        beta:       95
         |S|:   2^72.0
        prob: 0.125126
       scale:   12.760
          pp:       11
           d:     1666
      repeat:       35
         tau:      234
```

# The Hybrid Dual Attack

The Hybrid Dual Attack [CHHS19] adapts Albrecht's variant of the Dual attack on small and sparse secret LWE [Alb17] via the addition of a meet-in-the-middle process. In [CHHS19], this attack is shown to work particularly well for sparse secrets. The script accompanying [CHHS19] enables the cost of the Hybrid Dual attack to be estimated in the case of sparse secrets, but we found it too slow to be used for uniformly-random ternary secrets. 

To help solve this problem, we provide conservative runtime estimations (i.e., explicit underestimates) by assuming a square-root speed-up in the guessing phase of Albrecht's attack. This simulates the addition of a meet-in-the-middle process, under the same assumptions as the Hybrid Decoding attack above. To produce these estimates we use an adapted variant of the lwe-estimator [APS15,Est20]. Currently, this estimate only triggers for ternary secrets.

We briefly describe intuition for this attack. We begin by splitting the LWE matrix <img src = "images/lwe_matrix.svg" class = "inline"> into <img src = "images/lwe_matrix_split.svg" class = "inline"> where <img src = "images/A0_size.svg" class = "inline"> and <img src = "images/A1_size.svg" class = "inline">.

We perform lattice reduction on the basis of the lattice <img src = "images/lambda.svg" class = "inline">, where:

<img src = "images/lambda_description.svg" class = "inline">

Having retrieved a short vector <img src = "images/v_in_lambda.svg" class = "inline">, we check the inner product

<img src = "images/inner_product.svg" class = "inline">

and, if our initial guess of <img src = "images/tau.svg" class = "inline"> zeros is correct, then we can distinguish LWE from random using a collection of inner products as above. Otherwise, we can perform a guessing strategy over the dimension <img src = "images/tau.svg" class = "inline"> guessing space to correct for mis-guessed zeros. This technique is referred to as postprocessing in the LWE-Estimator [Est20]. In [Alb17, CHHS19] further details are given as to how to amortise the cost of lattice reduction (allowing for many "slightly longer" short vectors to be produced at a lower cost), as well as utilising scaled-normal form. 

Thus, the Hybrid Dual attack contains three phases:

1. Fixing the guessing dimension <img src = "images/tau.svg" class = "inline"> and constructing a basis for the lattice <img src = "images/lambda.svg" class = "inline"> defined above.
2. Performing lattice reduction (i.e. BKZ) on the basis of <img src = "images/lambda.svg" class = "inline">.
3. Performing batches of inner products (one batch for each guess), and checking the resulting distributions.

This guessing phase is typically realised using a meet-in-the-middle process, which requires exponential memory, but allows for the searching phase to be performed in a more (time) efficient manner. To simulate this meet-in-the-middle process, we conservatively assume a square-root speed-up, and ignore any probability loss introduced. Thus, the estimates provided by our script are explicit underestimates of security, based on the assumptions outlined below. 

## Assumptions

| Assumption | Choice | Notes |
|---|---|---|
| Reduced basis shape | Geometric Series Assumption |   |
| Meet in the middle speed-up | square-root |  |
| Meet in the middle success probability | 1 |  |
| Whether to consider cost of memory | No |  |

## What the Code Does

The code optimizes the complexity of the Hybrid Dual attack over the parameters <img src = "images/tau.svg" class = "inline"> and <img src = "images/beta.svg" class = "inline">, in the same manner that the LWE-Estimator [Est20] optimises Albrecht's dual attack, see [drop_and_solve](https://bitbucket.org/malb/lwe-estimator/src/a27675558616c19318e8df217ba4abedfe2fce07/estimator.py#lines-1677) and [dual_scale](https://bitbucket.org/malb/lwe-estimator/src/a27675558616c19318e8df217ba4abedfe2fce07/estimator.py#lines-2475). To retrieve estimates for the Hybrid Dual attack, we use a modified version of the LWE Estimator. There are very few changes, since the only addition we make is a square-root speed-up in the search phase of Albrecht's attack. These changes can be found [here](https://bitbucket.org/bencrts/lwe-estimator/commits/8f1e9b7149ec88ef2686b87aa22a754fbb493eba).

## Example usage of the code


```
%attach example_params.py
from hybrid_dual import estimator

estimate_lwe(reduction_cost_model = BKZ.sieve, **example_64())

         rop:   2^59.0
           m:      782
         red:   2^58.5
     delta_0: 1.009694
        beta:       91
      repeat:   2^13.6
           d:     1536
           c:   10.949
           k:      270
 postprocess:       16
```

We also provide an example which can be compared to the output of the [CHHS19] code. Their example is given at <url> https://github.com/swanhong/HybridLWEAttack#how-to-run</url>.

```
%attach example_params.py
from hybrid_dual import estimator

estimate_lwe(reduction_cost_model = BKZ.sieve, **chhs_19())

         rop:  2^120.5
           m:   2^11.7
         red:  2^120.0
     delta_0: 1.006151
        beta:      207
      repeat:   2^57.7
           d:   2^12.7
           c:   23.503
           k:   2^12.2
 postprocess:       15
```

This example is illustrative, and we note that the assumptions considered in [CHHS19] are vastly different for those considered here: namely, we do not consider probability loss in the meet-in-the-middle phase and we do not bound memory usage.


# Alternative Code and References

## Alternative code

There are several other implementations of hybrid attack estimates, we note these below. Feel free to add any other sources which we have missed.

### Hybrid Decoding

1. The NTRUPrime team's scripts: <url> https://ntruprime.cr.yp.to/estimate-20190329.sage </url>
2. Thomas Wunderer's scripts, hosted by Leo Ducas: <url>https://github.com/lducas/LatRedHybrid </url>
3. Rachel Player's edited variant of Wunderer's scripts: <url> https://github.com/rachelplayer/LatRedHybrid </url>
4. The script on Yongha Son's github page <url> https://github.com/Yongyongha/SparseLWE-estimator </url> associated to [SC19]

### Hybrid Dual

1. The script provided by [CHHS19]: <url>https://github.com/swanhong/HybridLWEAttack</url>
2. The script on Yongha Son's github page <url> https://github.com/Yongyongha/SparseLWE-estimator </url> associated to [CHHS19]

## References

[ACW19] <url>https://eprint.iacr.org/2019/1122</url> <br>
[APS15] <url>https://eprint.iacr.org/2015/046</url> <br>
[CHHS19] <url>https://eprint.iacr.org/2019/1114</url> <br>
[CP19] <url>https://eprint.iacr.org/2019/1148</url> <br>
[Est20] <url>https://bitbucket.org/malb/lwe-estimator</url><br>
[How07] <url>https://www.iacr.org/archive/crypto2007/46220150/46220150.pdf </url><br>
[SC19] <url> https://eprint.iacr.org/2019/1019.pdf </url>


