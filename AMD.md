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

## Code coverage
To produce a code coverage report for `ksw2_test`, first compile with coverage enabled:
```
make coverage=y
```
Run the executable one or more times
```
./ksw2_test .. 
```
Each time executable runs, code coverage data will be collected. To compile this data into a html
report, use:
```
make report
```
This report is written to the directory `coverage`.

## Code coverage script 
The script `coverage.sh` is designed to run `ksw2_test` multiple times to take all branches inside a
kernel. This
script can be used in conjunction with the `run_test_suite` script as follows:
```
 KSW2ID=extd2_cpp KSW2=./coverage.sh ./run_test_suite.sh test_suite -t extd2_cpp # test
 KSW2ID=extd2_sse KSW2=./coverage.sh ./run_test_suite.sh test_suite -t extd2_sse # reference
 ```
 Set `KSW2COV` if the main executable is not `./ksw2-test`.
 Use the comparison script as before,
 ```
 ./compare_test_suite.sh extd2_cpp extd2_sse
 ```
 To design your own code coverage tests, put one or more calls to `ksw2-test` in a script. See
 `coverage.sh` for an example.



