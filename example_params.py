import sys
sys.path.insert(0,'../estimator')
import estimator as est

## Example parameter sets

def example_64():
    q = 2**47
    return {"n": 1024,
            "alpha": est.alphaf(3.19, q, True),
            "q": q,
            "m": 1024,
            "secret_distribution": ((-1, 1), 64)}

def example_binary_64():
    q = 2**47
    return {"n": 1024,
            "alpha": est.alphaf(3.19, q, True),
            "q": q,
            "m": 1024,
            "secret_distribution": (0, 1)}

def example_128():
    q = 2**47
    return {"n": 1024,
            "alpha": est.alphaf(3.19, q, True),
            "q": q,
            "m": 1024,
            "secret_distribution": ((-1, 1), 128)}

def example_ternary():
    q = 2**200
    return {"n": 4096,
            "alpha": est.alphaf(3.19, q, True),
            "q": q,
            "m": 1024,
            "secret_distribution": (-1, 1)}

def chhs_19_repository_example():
    q = 2**125
    return {"n": 8192,
            "alpha": est.alphaf(3.19, q, True),
            "q": q,
            "m": 8192,
            "secret_distribution": ((-1, 1),64)}

def ntruprime():
    q = 4591
    return {"n": 761,
            "alpha": est.alphaf(sqrt(2./3), q, True),
            "q": q,
            "m": 761,
            "secret_distribution": ((-1, 1),250)}

def tfhe():
    q = 2**32
    return {"n": 1024,
            "alpha": sqrt(2*pi)*2**(-25),
            "q": q,
            "m": 1024,
            "secret_distribution": (0,1)}
