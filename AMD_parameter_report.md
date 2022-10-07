# Code coverage

This report discusses some of the modifications that we (AMD) have done to KSW2 test in order to
set additional flags and parameters by passing input arguments. We also discuss all of the flags and
parameters that we have tested.

## Flags

### KSW_EZ_SCORE_ONLY

`KSW_EZ_SCORE_ONLY` is not used in minimap2 but it is used in ksw2. This option disables the use of
"cigar". Currently, the `extd2_cpp` kernel segmentation faults.

```
./ksw2-test -t extd2_cpp test_suite/small/q0.fa test_suite/small/t0.fa -s
```

GDB debug log:
```
Program received signal SIGSEGV, Segmentation fault.
0x0000555555570d6c in ksw_extd2_cpp (km=0x0, qlen=45, query=0x5555555970e0 "", tlen=45, target=0x555555597140 "", m=5 '\005',
    mat=0x7fffffffdc10 "\002\374\374", <incomplete sequence \374>, q=4 '\004', e=2 '\002', q2=13 '\r', e2=1 '\001', w=45, zdrop=-1, end_bonus=0,
    flag=1, ez=0x7fffffffdbd0) at ksw2_extd2.cpp:533
533             off[r] = st;
```

### KSW_EZ_RIGHT
This flag is used in minimap and is used for right-alignment. It is enabled by passing `-r` to
kws2_test.

### KSW_EZ_GENERIC_SC
This flag is not used in minimap2. It is also not available as option to ksw2_test. ksw2_test has been
modified to test this flag. Pass `-c` to enable it.
```
./ksw2-test -t extd2_cpp test_suite/small/q0.fa test_suite/small/t0.fa -c
```

### KSW_EZ_APPROX_MAX and KSW_EZ_APPROX_DROP
These flags are not supported by the GPU implementation
* KSW_EZ_APPROX_MAX: Approximate max operation to improve performance (performance optimization). 
* KSW_EZ_APPROX_DROP: Approximate Z-drop to improve performance (performance optimization).

In ksw2_test, both of these flags are enabled by passing `-g`. 
```
./ksw2-test -t extd2_sse test_suite/small/q0.fa test_suite/small/t0.fa -g
```

### KSW_EZ_EXTZ_ONLY
Use extension only (backtracking). This option is not available in ksw2, ksw2_test has been
modified to test this flag. Pass -x` to enable it.

```
./ksw2-test -t extd2_cpp test_suite/small/q0.fa test_suite/small/t0.fa -x
```

## Parameters

* end_bonus : This parameter is used during the backtracking stage and can be used to trigger the
  code path:
```
		} else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
		ez->reach_end = 1;
		ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
```
* -E : Increases the gap penalty can cause `long_thres` to be incremented.
