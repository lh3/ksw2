import numpy as np

np.random.seed(0)

# Convert integer to base
def tobase(x):
    bases = ['G', 'C', 'T', 'A', 'g', 'c', 't', 'a']
    assert(x < 8)
    return bases[x]


# Convert base to integer
def toint(x):
    bases = {'G' : 0, 'C' : 1, 'T' : 2, 'A' : 3, 'g' : 4, 'c': 5, 't' : 6, 'a': 7} 
    return bases[x]


def gen_random_query(num_queries=1, min_len=1, max_len=10000, verbose=False):
    """
    write_random_query writes a ksw2 input file with randomized strings.

    Args:
        filename: The name of the output file
        num_queries: The number of queries (maps to rows) in the output
        min_len: The minimum length of each query string
        max_len: The maximum length of each query string
        verbose: Give verbose output.

    """

    if verbose:
        print(f"Settings num_queries={num_queries}, min_len={min_len}, max_len={max_len}")

    alphabet = np.ones(max_len) * 8
    q = []

    assert(min_len <= max_len)

    for i in range(num_queries):
        ql = np.random.randint(max_len - min_len + 1) + min_len
        assert(ql >= min_len)
        assert(ql <= max_len)
        qi = np.random.randint(alphabet[:ql])
        seq = map(tobase, qi)
        qbase = np.array(list(seq))
        q.append(qbase)

    return q


def perturb_query(q, threshold):
    """
    Perturb each entry in a list of queries based on a given `threshold`.
    Only entries that by chance exceed the chosen threshold are modified.


    Args:
        q: List of queries
        threshold: The probability that an entry is modified (0 <= t <= 1)

    Returns:
        dq: List of perturbed queries

    """
    dq = []

    assert(threshold <= 1.0)
    assert(threshold >= 0.0)

    max_len = np.max([len(qi) for qi in q])

    alphabet = np.ones(max_len) * 8

    for qi in q:
        # Convert to integer representation
        seq = map(toint, qi)
        qint = np.array(list(seq))

        # Modify the entries in qint that falls below the threshold
        p = np.random.random(len(qint))
        pi = np.random.randint(alphabet[:len(qint)])
        mask = p < threshold
        qint[mask] = pi[mask]

        # Convert back to base presentation
        qbase = np.array(list(map(tobase, qint)))

        dq.append(qbase)


    return dq

def cut_query(q, max_cut):
    """
    Make deletions to a list of queries randomly. 
    
    Args:
        q: List of queries
        max_cut: [ max cuts to make, max length of each cut w.r.t. query length  ( 0 <= lmac < 1 )]

    Returns:
        cq: List of cut queries.
    """
    cq = []
    for qi in q:
        c_count = np.random.randint(max_cut[0])
        for i in range(c_count):
            cl = round(len(qi) * max_cut[1])
            cl = np.random.randint(cl+1)
            st = np.random.randint(len(qi) - cl)
            en = st + cl
            assert(en < len(qi))
            
            qi = np.concatenate((qi[:st], qi[en:]))
        cq.append(qi)
    return cq


def write_ksw2file(filename, label, queries, verbose=False):
    """
    write_ksw2file writes a ksw2 input file.

    Args:
        filename: The name of the output file
        label: The label of the query variable, typically "q" or "t"
        queries: A list of queries (can be generated using gen_random_query)

    """


    num_queries = len(queries)

    with open(filename, "w") as fh:
        
        for i, query in enumerate(queries):
            num = i + 1
            fh.write(f">{label}{num}\n")
            fh.write("".join(query) + "\n")

    if verbose:
        print("Wrote:", filename)


def write_test_case(path, num_queries=10, min_len=1, max_len=100, thresholds=[0.1], max_cut=[2, 0.1], verbose=False):
    import os
    try:
        os.makedirs(path)
    except:
        pass

    if verbose:
        print(f"Generating test case {path} ...")

    for num, p in enumerate(thresholds):
        q = gen_random_query(num_queries, min_len, max_len, verbose)
        t = perturb_query(q, p)
        q = cut_query(q, max_cut)
        t = cut_query(t, max_cut)

        write_ksw2file(f"{path}/q{num}.fa", "q", q, verbose)
        write_ksw2file(f"{path}/t{num}.fa", "t", t, verbose)
