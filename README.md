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

* [ksw2_gg.c][ksw2_gg.c]: global alignment; Green's standard formulation
* [ksw2_gg2.c][ksw2_gg2.c]: global alignment; Suzuki's diagonal formulation
* [ksw2_gg2_sse_u.c][ksw2_gg2_sse_u.c]: global alignment with unaligned SSE intrinsics; Suzuki's
* [ksw2_gg2_sse.c][ksw2_gg2_sse.c]: global alignment with mostly aligned SSE intrinsics; Suzuki's
* [ksw2_extz.c][ksw2_extz.c]: alignment extension; Green's formulation
* [ksw2_extz2_sse_u.c][ksw2_extz2_sse_u.c]: extension with unaligned SSE intrinsics; Suzuki's
* [ksw2_extz2_sse.c][ksw2_extz2_sse.c]: extension with mostly aligned SSE intrinsics; Suzuki's

Users are encouraged to copy the header file `ksw2.h` and relevant
`ksw2_*.c` file to their own source code trees. On x86 CPUs with SSE2
intrinsics, `ksw2_extz2_sse.c` is recommended in general. It supports global
alignment, alignment extension with Z-dropoff, score-only alignment and
right-aligned CIGARs. `ksw2_gg2_sse.c` is faster than `ksw2_extz2_sse.c` for
pure left-aligned global alignment.

## Limitations

* No dynamic banding. [edlib][edlib] dynamically changes the band width to
  always find the optimal alignment without traversing every cell in the DP
  matrix. The same strategy can be applied to KSW2 in principle. However, in
  case of alignment extension, dynamic banding might not help much.

* No adaptive banding. [libgaba][gaba] implements another type of dynamic
  banding called [adaptive banding][adap-band] (which I have thought of
  independently). It moves the band of fixed width based on the score at the
  current anti-diagonal. It effectively increases the band width at little
  computational cost. However, this strategy does not help much to find
  long gaps. It is not clear about its added value without more careful
  evaluations.

* The inner loop with SSE intrinsics could probably be improved. On a few
  examples, `ksw2_gg2` is about twice as slow as `ksw2_gg`, suggesting the
  diagonal formulation is less efficient. At the same time, `ksw2_extz2_sse` is
  four times as fast as `ksw2_extz`, suggesting SSE works, but the speedup
  is not as good as 8 or even 16 fold. There may be room for improvement.

* No thorough test or benchmarking as of now.

* No AVX implementations.

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
