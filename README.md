## What

KSW2 is a library to align a pair of biological sequences based on dynamic
programming (DP). So far it comes with global alignment and alignment extension
(no local alignment yet) under affine gap penalty. It supports fixed banding
and optionally produces alignment paths (i.e. CIGARs) with gaps either left- or
right-aligned.  In addition to plain implementations of the algorithms, KSW2
also provides implementations using SSE2 and SSE4.1 intrinsics. It adopted
[Hajime Suzuki][hs]'s [formulation][hs-eq] which enables 16-way SSE
parallelization for the most part of the inner loop, regardless of the maximum
score of the alignment.

## Why

Many libraries perform DP-based alignment. The following table gives an
overview:

|Library         |CIGAR|Intra-seq|Affine-gap|Local    |Global   |Glocal   |Extension|
|:---------------|:---:|:-------:|:--------:|:-------:|:-------:|:-------:|:-------:|
|[edlib][edlib]  |Yes  |Yes      |No        |Very fast|Very fast|Very fast|N/A      |
|[KSW][klib]     |Yes  |Yes      |Yes       |Fast     |Slow     |N/A      |Slow     |
|KSW2            |Yes  |Yes      |Yes       |N/A      |Fast     |N/A      |Fast     |
|[libgaba][gaba] |Yes  |Yes      |Yes       |N/A?     |N/A?     |N/A?     |Fast     |
|[libssa][ssa]   |No   |No?      |Yes       |Fast     |Fast     |N/A      |N/A      |
|[Opal][opal]    |No   |No       |Yes       |Fast     |Fast     |Fast     |N/A      |
|[Parasail][para]|No   |Yes      |Yes       |Fast     |Fast     |Fast     |N/A      |
|[SeqAn][seqan]  |Yes  |Yes      |Yes       |Slow     |Slow     |Slow     |N/A      |
|[SSW][ssw]      |Yes  |Yes      |Yes       |Fast     |N/A      |N/A      |N/A      |
|[SWIPE][swipe]  |Yes  |No       |Yes       |Fast     |N/A?     |N/A?     |N/A      |
|[SWPS3][swps3]  |No   |Yes      |Yes       |Fast     |N/A?     |N/A      |N/A      |

We developed KSW2 because it comes with a set of features needed for developing
aligners. For a seed-and-extend based aligner, KSW2 can be used to close gaps
between seeds.  When there is a long gap between two adjacent seed hits, we
prefer a global alignment but cannot force the query and reference sequences to
be fully aligned because they may differ due to structural variations such as long
inversions, which may leave a long poorly aligned regions in the middle.  KSW2
can detect such poorly regions with *diagonal X-drop*, which I call as *Z-drop*.
Z-drop is like X-drop except that it does not panelize gap extensions and thus
helps to recover long gaps. Variant callers for high-throughput sequencing data
usually expect gaps to be left-aligned.  To achieve this, we need gaps to be
left-aligned when we extend to the right, while right-aligned when we extend to
the left. Both gap placements are necessary.

## How to use

Each `ksw2_*.c` file implements a single function and is independent of each
other. Here are brief descriptions about what each file implements:

* [ksw2_gg.c](ksw2_gg.c): global alignment; Green's standard formulation
* [ksw2_gg2.c](ksw2_gg2.c): global alignment; Suzuki's diagonal formulation
* [ksw2_gg2_sse.c](ksw2_gg2_sse.c): global alignment with mostly aligned SSE intrinsics; Suzuki's
* [ksw2_extz.c](ksw2_extz.c): alignment extension; Green's formulation
* [ksw2_extz2_sse.c](ksw2_extz2_sse.c): extension with mostly aligned SSE intrinsics; Suzuki's

Users are encouraged to copy the header file `ksw2.h` and relevant
`ksw2_*.c` file to their own source code trees. On x86 CPUs with SSE2
intrinsics, `ksw2_extz2_sse.c` is recommended in general. It supports global
alignment, alignment extension with Z-drop, score-only alignment, global-only
alignment and right-aligned CIGARs. `ksw2_gg*.c` are mostly for demonstration
and comparison purposes. They are annotated with more comments and easier to
understand than `ksw2_ext*.c`. Header file [ksw2.h](ksw2.h) gives brief
documentations.

To compile the test program `ksw-test`, just type `make`. For most x86
compilers, this doesn't take the advantage of SSE4.1. To compile with SSE4.1
for better performance, use `make sse4=1` instead. If you have installed
parasail, use `make sse4=1 parasail=prefix`, where `prefix` points to the
parasail install directory (e.g. `/usr/local`). The test program can also
use libgaba, but I failed to produce desired alignment on the test data,
probably due to my fault.

## Performance Analysis

The following table shows timing on two pairs of long sequences (both in the
"test" directory).

|Data set|Command line options             |Time (s)|CIGAR|Ext|SIMD|Source  |
|:-------|:--------------------------------|:-------|:---:|:-:|:--:|:-------|
|50k     |-t gg -s                         |7.3     |N    |N  |N   |ksw2    |
|        |-t gg2 -s                        |19.8    |N    |N  |N   |ksw2    |
|        |-t extz -s                       |9.2     |N    |Y  |N   |ksw2    |
|        |-t ps\_nw                        |9.8     |N    |N  |N   |parasail|
|        |-t ps\_nw\_striped\_sse2\_128\_32|2.9     |N    |N  |SSE2|parasail|
|        |-t ps\_nw\_striped\_32           |2.2     |N    |N  |SSE4|parasail|
|        |-t ps\_nw\_diag\_32              |3.0     |N    |N  |SSE4|parasail|
|        |-t ps\_nw\_scan\_32              |3.0     |N    |N  |SSE4|parasail|
|        |-t extz2\_sse -sg                |0.96    |N    |N  |SSE2|ksw2    |
|        |-t extz2\_sse -sg                |0.84    |N    |N  |SSE4|ksw2    |
|        |-t extz2\_sse -s                 |3.0     |N    |Y  |SSE2|ksw2    |
|        |-t extz2\_sse -s                 |2.7     |N    |Y  |SSE4|ksw2    |
|16.5k   |-t gg -s                         |0.84    |N    |N  |N   |ksw2    |
|        |-t gg                            |1.6     |Y    |N  |N   |ksw2    |
|        |-t gg2                           |3.3     |Y    |N  |N   |ksw2    |
|        |-t extz                          |2.0     |Y    |Y  |N   |ksw2    |
|        |-t extz2\_sse                    |0.40    |Y    |Y  |SSE4|ksw2    |
|        |-t extz2\_sse -g                 |0.18    |Y    |N  |SSE4|ksw2    |

The standard DP formulation is about twice as fast as Suzuki's diagonal
formulation (`-tgg` vs `-tgg2`), but SSE-based diagonal formulation
is several times faster than the standard DP. If we only want to compute one
global alignment score, we can use 16-way parallelization in the entire inner
loop.  For extension alignment, though, we need to keep an array of 32-bit
scores and have to use 4-way parallelization for part of the inner loop. This
significantly reduces performance (`-sg` vs `-s`).  KSW2 is faster than
parasail partly because the former uses one score for all matches and another
score for all mismatches. For diagonal formulations, vectorization is more
complex given a generic scoring matrix.

It is possible to further accelerate global alignment with dynamic banding as
is implemented in [edlib][edlib]. However, it is not as effective for extension
alignment. Another idea is [adaptive banding][adap-band], which might be worth
trying at some point.



[hs]: https://github.com/ocxtal
[hs-eq]: https://github.com/ocxtal/diffbench
[edlib]: https://github.com/Martinsos/edlib
[klib]: https://github.com/attractivechaos/klib
[para]: https://github.com/jeffdaily/parasail
[opal]: https://github.com/Martinsos/opal
[ssw]: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
[ssa]: https://github.com/RonnySoak/libssa
[gaba]: https://github.com/ocxtal/libgaba
[adap-band]: https://github.com/ocxtal/adaptivebandbench
[swipe]: https://github.com/torognes/swipe
[swps3]: http://lab.dessimoz.org/swps3/
[seqan]: http://seqan.de
