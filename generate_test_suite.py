import test_suite as test

test.write_test_case("test_suite/small", num_queries=100, min_len=1, max_len=100, thresholds=[0.1, 0.5,
    0.75, 0.9], verbose=True)

test.write_test_case("test_suite/medium", num_queries=100, min_len=100, max_len=1000, thresholds=[0.1, 0.5,
    0.75, 0.9], verbose=True)

test.write_test_case("test_suite/large", num_queries=100, min_len=1000, max_len=20000, thresholds=[0.1, 0.5,
    0.75, 0.9], verbose=True)

#test.write_test_case("test_suite/huge", num_queries=4, min_len=50000, max_len=100000, thresholds=[0.1],
#                     verbose=True)
#
#test.write_test_case("test_suite/giant", num_queries=4, min_len=100000, max_len=4000000, thresholds=[0.1],
#                     verbose=True)
