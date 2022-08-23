# KSW2 Test suite

1. Clone and modify `cli.c` from https://github.com/lh3/ksw2 so that you can call your own kernels,
   e.g., `extd2`.
Compiling this code produces the executable `ksw2-test`.

2. Generate a test suite using
```bash
python3 generate_test_suite.py
```
This script will generate random data with a fixed seed for different alignment lengths. Modify this
script to build your own test suites. In the current test suite, the "giant" test crashes
`ksw2-test` due to insufficient memory. To overcome this issue, specify the bandwidth `-w` argument. 

3. Once you have generated data, run:
```bash
./run_test_suite.sh
```
Simply running the command in bash will show usage instructions.For example,
```bash
KSW2ID=reference ./run_test_suite.sh test_suite -t extd2_sse
```
will run `ksw2_test` for all tests in `test_suite` using the kernel `extd2_sse` (this name depends
on how `cli.c` is configured). Run it using your own kernel by passing the name of your kernel to
run script, e.g.,

```bash
./run_test_suite.sh test_suite -t extd2
```
By omitting, `KSW2ID`, the script will automatically create a output directory for you.

4. Given both reference runs and your own runs, compare the two via
```bash
./compare_test_suite.sh myrun reference
```
