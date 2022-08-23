import numpy as np
import test_suite as test

def test_perturb_query():
    q = test.gen_random_query(max_len=10)
    p = test.perturb_query(q, 0.0)
    for qi, pi in zip(q, p):
        assert np.all(qi == pi)

    err = 0.05
    t = 0.1
    q = test.gen_random_query(min_len=10000, max_len=10000)
    p = test.perturb_query(q, t)

    for qi, pi in zip(q, p):
        assert abs(np.sum(qi == pi) / len(qi) - (1 - t) ) < err


