## What

KSW2 is a library to align a pair of biological sequences based on dynamic
programming (DP). So far it comes with global alignment and alignment extension
(no local alignment yet) under affine gap penalty. It supports fixed banding
and optionally produces alignment paths (i.e. CIGARs) with gaps either left- or
right-aligned.  In addition to plain implementations of the algorithms, KSW2
also provides implementations using SSE2 and SSE4.1 intrinsics. It adopted
[Hajime Suzuki][hs]'s [formulation][hs-eq] which enables 16-way SSE
parallelization regardless of the maximum score of the alignment.

## Why

KSW2 comes with a set of features needed for developing aligners. For a
seed-and-extend based aligner, KSW2 can be used to close gaps between seeds.
When there is a long gap between two adjacent seed hits, we prefer a global
alignment but cannot force the query and reference sequences to be aligned
because they may differ due to structural variations such as long inversions.
KSW2 can detect poorly aligned regions with diagonal X-drop, which I call as
Z-drop. Z-drop is like X-drop except that it does not panelize gap extensions
and thus helps to recover long gaps. Variant callers for high-throughput
sequencing data usually expect gaps to be left-aligned.  To achieve this, we
need gaps to be left-aligned when we extend to the right, while right-aligned
when we extend to the left.

Many libraries perform DP-based alignment. [KSW][klib], the predecessor of
KSW2, implements the same set of features, but it has minor bugs and is slow
for global alignment and alignment extension.  [edlib][edlib] is very fast but
it does not support affine gap penalty and requires the alignment to reach the
end of the query. [Parasail][para] and [libssa][ssa] do not produce CIGARs or
support extension. [Opal][opal] implements inter-sequence parallelization which
is not applicable to my use cases. [SSW][ssw] comes with local alignment only.
[libgaba][gaba] is probably the closest to KSW2, but it still lacks a few
subtle features such as diagonal X-drop.

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
understand than `ksw2_ext*.c`.

## Performance Analysis

The following table shows timing on two pairs of sequences (both in the "test"
directory).

|Data set|Command line options             |Time (s)|CIGAR|Ext|SIMD|Source  |
|:-------|:--------------------------------|:-------|:---:|:-:|:--:|:-------|
|50k     |-t gg -s                         |7.3     |N    |N  |N   |ksw2    |
|        |-t extz -s                       |9.2     |N    |Y  |N   |ksw2    |
|        |-t ps\_nw                        |9.8     |N    |N  |N   |parasail|
|        |-t ps\_nw\_striped\_sse2\_128\_32|2.9     |N    |N  |SSE2|parasail|
|        |-t ps\_nw\_striped\_32           |2.2     |N    |N  |SSE4|parasail|
|        |-t ps\_nw\_diag\_32              |3.0     |N    |N  |SSE4|parasail|
|        |-t ps\_nw\_scan\_32              |3.0     |N    |N  |SSE4|parasail|
|        |-t extz2\_sse -s                 |3.0     |N    |Y  |SSE2|ksw2    |
|        |-t extz2\_sse -s                 |2.7     |N    |Y  |SSE4|ksw2    |
|        |-t extz2\_sse -sg                |0.96    |N    |N  |SSE2|ksw2    |
|        |-t extz2\_sse -sg                |0.84    |N    |N  |SSE4|ksw2    |
|16.5k   |-t gg -s                         |0.84    |N    |N  |N   |ksw2    |
|        |-t gg                            |1.6     |Y    |N  |N   |ksw2    |
|        |-t gg2                           |3.3     |Y    |N  |N   |ksw2    |
|        |-t extz                          |2.0     |Y    |Y  |N   |ksw2    |
|        |-t extz2\_sse                    |0.40    |Y    |Y  |SSE4|ksw2    |
|        |-t extz2\_sse -g                 |0.18    |Y    |N  |SSE4|ksw2    |

The standard DP formulation is about twice as fast as Suzuki's diagonal
formulation (`-tgg` vs `-tgg2`), but SSE4-based diagonal formulation
is several times faster than the standard DP. If we only want to compute one
global alignment score, we can use 16-way parallelization throughout.  For
extension alignment, though, we need to keep an array of 32-bit scores, which
significantly reduces performance (`-sg` vs `-s`).  KSW2 is faster than
parasail partly because the former uses one score for all matches and another
score for all mismatches. Vectorization is harder given a generic scoring
matrix.

It is possible to further accelerate global alignment with dynamic banding as
is implemented in [edlib][edlib]. However, it might not be as effective for
extension alignment. Another idea is [adaptive banding][adap-band], which
might be worth trying at some point.



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
