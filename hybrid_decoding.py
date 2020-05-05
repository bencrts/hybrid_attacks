# -*- coding: utf-8 -*-
"""
Seucrity estimates for the Hybrid Decoding attack

Requires a local copy of the LWE Estimator from https://bitbucket.org/malb/lwe-estimator/src/master/estimator.py in a folder called "estimator"
"""

from sage.all import ZZ, binomial, sqrt, log, exp, oo, pi, prod, RR
from sage.probability.probability_distribution import RealDistribution
from estimator import estimator as est

## Core cost models

core_sieve = lambda beta, d, B: ZZ(2)**RR(0.292*beta + 16.4)
core_qsieve = lambda beta, d, B: ZZ(2)**RR(0.265*beta + 16.4)

## Utility functions

def _primal_scale_factor(secret_distribution, alpha=None, q=None, n=None):
    """
    Scale factor for primal attack, in the style of [BaiGal14]. In the case of non-centered secret
    distributions, it first appropriately rebalances the lattice bases to maximise the scaling.
    Note that this function has been copied from the LWE-Estimator as importing threw an error.

    :param secret_distribution: distribution of secret, see module level documentation for details
    :param alpha: noise rate `0 ≤ α < 1`, noise has standard deviation `αq/sqrt{2π}`
    :param q: modulus `0 < q`
    :param n: only used for sparse secrets

    """

    # For small/sparse secret use Bai and Galbraith's scaled embedding
    # NOTE: We assume a <= 0 <= b

    scale = RR(1)
    if est.SDis.is_small(secret_distribution):
        # target same same shortest vector length as in Kannan's embedding of same dimension
        stddev = est.stddevf(alpha*q)
        var_s = est.SDis.variance(secret_distribution, alpha, q, n=n)
        if stddev**2 > var_s:
            # only scale if the error is sampled wider than the secret
            scale = stddev/RR(sqrt(var_s))

    return scale


def max_blocksize(secbits, reduction_cost_model = est.BKZ.sieve):
    """ A function to determine the maximal blocksize using CORE cost models.
        This provides a (strict) upper bound on beta for our computations, to allow
        for more efficient searching
    :param secbits: the target security level
    :param reduction_cost_model: the lattice reduction cost model used for BKZ

    """

    rop = 0
    beta_max = 40
    while rop <= secbits:
        beta_max += 5
        rop = log(est.lattice_reduction_cost(reduction_cost_model, est.delta_0f(beta_max), 1)["rop"],2)

    return beta_max


def sq_GSO(d, beta, det):
    """
    Return squared GSO lengths after lattice reduction according to the GSA

    :param q: LWE modulus
    :param d: lattice dimension
    :param beta: blocksize used in BKZ
    :param det: lattice determinant

    """

    r = []
    for i in range(d):
        r_i = est.delta_0f(beta)**(((-2*d*i) / (d-1))+d) * det**(1/d)
        r.append(r_i**2)

    return r


def babai_probability_wun16(r, norm):
     """
     Compute the probability of Babai's Nearest Plane, using techniques from the NTRULPrime submission to NIST

     :param r: squared GSO lengths
     :param norm: expected norm of the target vector

     """
     R = [RR(sqrt(t)/(2*norm)) for t in r]
     T = RealDistribution('beta', ((len(r)-1)/2,1./2))
     probs = [1 - T.cum_distribution_function(1 - s**2) for s in R]
     return prod(probs)


## Estimate hybrid decoding attack complexity

def hybrid_decoding_attack(n, alpha, q, m, secret_distribution,
                   beta, tau = None, mitm=True, reduction_cost_model=est.BKZ.sieve):
    """
    Estimate cost of the Hybrid Attack,

    :param n: LWE dimension `n > 0`
    :param alpha: noise rate `0 ≤ α < 1`, noise will have standard deviation `αq/sqrt{2π}`
    :param q: modulus `0 < q`
    :param m: number of LWE samples `m > 0`
    :param secret_distribution: distribution of secret
    :param beta: BKZ block size β
    :param tau: guessing dimension τ
    :param mitm: simulate MITM approach (√ of search space)
    :param reduction_cost_model: BKZ reduction cost model

    EXAMPLE:

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

    """

    n, alpha, q = est.Param.preprocess(n, alpha, q)

    # d is the dimension of the attack lattice
    d = m + n - tau

    # h is the Hamming weight of the secret
    # NOTE: binary secrets are assumed to have Hamming weight ~n/2, ternary secrets ~2n/3
    # this aligns with the assumptions made in the LWE Estimator
    h = est.SDis.nonzero(secret_distribution, n=n)
    sd = alpha*q/sqrt(2*pi)

    # compute the scaling factor used in the primal lattice to balance the secret and error
    scale = _primal_scale_factor(secret_distribution, alpha=alpha, q=q, n=n)

    # 1. get squared-GSO lengths via the Geometric Series Assumption
    # we could also consider using the BKZ simulator, using the GSA is conservative
    r = sq_GSO(d, beta, q**m * scale**(n-tau))

    # 2. Costs
    bkz_cost = est.lattice_reduction_cost(reduction_cost_model, est.delta_0f(beta), d)
    enm_cost = est.Cost()
    enm_cost["rop"] = d**2/(2**1.06)

    # 3. Size of search space
    # We need to do one BDD call at least
    search_space, prob, hw = ZZ(1), 1.0, 0

    # if mitm is True, sqrt speedup in the guessing phase. This allows us to square the size
    # of the search space at no extra cost.
    # NOTE: we conservatively assume that this mitm process succeeds with probability 1.
    ssf = sqrt if mitm else lambda x: x

    # use the secret distribution bounds to determine the size of the search space
    a, b = est.SDis.bounds(secret_distribution)

    # perform "searching". This part of the code balances the enm_cost with the cost of lattice
    # reduction, where enm_cost is the total cost of calling Babai's algorithm on each vector in
    # the search space.

    if tau:
        prob = est.success_probability_drop(n, h, tau)
        hw = 1
        while hw < h and hw < tau:
            prob += est.success_probability_drop(n, h, tau, fail=hw)
            search_space += binomial(tau, hw) * (b-a)**hw

            if enm_cost.repeat(ssf(search_space))["rop"] > bkz_cost["rop"]:
                # we moved too far, so undo
                prob -= est.success_probability_drop(n, h, tau, fail=hw)
                search_space -= binomial(tau, hw) * (b-a)**hw
                hw -= 1
                break
            hw += 1

        enm_cost = enm_cost.repeat(ssf(search_space))

    # we use the expectation of the target norm. This could be longer, or shorter, for any given instance.
    target_norm = sqrt(m * sd**2 + h * RR((n-tau)/n) * scale**2)

    # account for the success probability of Babai's algorithm
    prob*=babai_probability_wun16(r, target_norm)

    # create a cost string, as in the LWE Estimator, to store the attack parameters and costs
    ret = est.Cost()
    ret["rop"] = bkz_cost["rop"] + enm_cost["rop"]
    ret["pre"] = bkz_cost["rop"]
    ret["enum"] = enm_cost["rop"]
    ret["beta"] = beta
    ret["|S|"] = search_space
    ret["prob"] = prob
    ret["scale"] = scale
    ret["pp"] = hw
    ret["d"] = d
    ret["tau"] = tau

    # 5. Repeat whole experiment ~1/prob times
    ret = ret.repeat(est.amplify(0.99, prob), select={"rop": True,
                                                      "pre": True,
                                                      "enum": True,
                                                      "beta": False,
                                                      "d": False,
                                                      "|S|": False,
                                                      "scale": False,
                                                      "prob": False,
                                                      "pp": False,
                                                      "tau": False})

    return ret


## Optimize attack parameters

def parameter_search(n, alpha, q, m, secret_distribution, mitm = True, reduction_cost_model=est.BKZ.sieve, secbits = None):

    """
    :param n: LWE dimension `n > 0`
    :param alpha: noise rate `0 ≤ α < 1`, noise will have standard deviation `αq/sqrt{2π}`
    :param q: modulus `0 < q`
    :param m: number of LWE samples `m > 0`
    :param secret_distribution: distribution of secret
    :param beta_search: tuple (β_min,  β_max, granularity) for the search space of β, default is (60,301,20)
    :param tau: tuple (τ_min, τ_max, granularity) for the search space of τ, default is (0,501,20)
    :param mitm: simulate MITM approach (√ of search space)
    :param reduction_cost_model: BKZ reduction cost model

    EXAMPLE:

    parameter_search(mitm = False, reduction_cost_model = est.BKZ.sieve, **example_64())

         rop:   2^69.5
         pre:   2^68.9
        enum:   2^68.0
        beta:      110
         |S|:   2^40.9
        prob: 0.045060
       scale:   12.760
          pp:        6
           d:     1730
      repeat:      100
         tau:      170

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

    """

    primald = est.partial(est.drop_and_solve, est.dual_scale, postprocess=True, decision=True)
    bl = primald(n, alpha, q, secret_distribution=secret_distribution, m=m, reduction_cost_model=reduction_cost_model)

    # we take the number of LWE samples used to be the same as in the primal attack in the LWE Estimator
    m = bl["m"]

    f = est.partial(hybrid_decoding_attack, n=n, alpha=alpha, q=q, m=m, secret_distribution=secret_distribution,
                    reduction_cost_model=reduction_cost_model,
                    mitm=mitm)

    best = None

    # when searching for parameters of a specific security level (e.g 128-bits), we can use the input parameter
    # 'secbits'. This gives us an upper bound on the blocksize being used (this will not search for any beta
    # such that the cost of BKZ(beta) on a dimension 1 lattice is > 2**secbits). This is useful when deriving
    # parameters for a given security level, but should not be used when finding the security level of a given
    # parameter set.

    if secbits is None:
        beta_max = n
    else:
        # compute max blocksize based on a dimension 1 lattice
        beta_max = max_blocksize(secbits, reduction_cost_model = reduction_cost_model)
        print("the maximal blocksize is {}".format(beta_max))

    print("BETA MAX IS {}".format(beta_max))

    # NOTE: we decribe our searching strategy below. To produce more accurate estimates,
    # change this part of the code to ensure a more granular search. As we are using
    # homomorphic-encryption style parameters, the running time of the code can be quite high,
    # justifying the below choices.
    # We start at beta = 60 and go up to beta_max in steps of 50

    beta_search = (60, beta_max, 50)

    for beta in range(beta_search[0], beta_search[1], beta_search[2])[::-1]:
        tau = 0
        best_beta = None
        while tau < n:
                cost = f(beta=beta, tau=tau)
                cost["tau"] = tau
                if best_beta is None:
                    best_beta = cost
                if RR(log(cost["rop"],2)) < RR(log(best_beta["rop"],2)):
                    best_beta = cost
                if best is None:
                    best = cost
                if RR(log(cost["rop"],2)) < RR(log(best["rop"],2)):
                    best = cost
                # tau is searched in steps of n//10
                tau += n//10

    # now do a second, more granular search
    # we start at the beta which produced the lowest running time, and search ± 25 in steps of 10
    for beta in range(best["beta"] - 25, best["beta"] + 25, 10)[::-1]:
        # we start at the tau which produced the lowest running time, and search ± n//20 in steps of n//100 or 1, whichever is higher
        for tau in range(max(best["tau"] - n//20,0), best["tau"] + n//20, max(n//100,1)):
            cost = f(beta=beta, tau=tau)
            cost["tau"] = tau
            if cost["rop"] < best["rop"]:
                best = cost
    return best
