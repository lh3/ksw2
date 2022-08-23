#!/usr/bin/bash
if [ -z "$KSW2" ]; then
    KSW2=./ksw2-test
fi

if [ -z "$KSW2ID" ]; then
    d=`date`
    hash=`echo -n $d | sha256sum | cut -c1-6`
    KSW2ID=$hash
fi

init() {
    echo "KSW2 Test Suite ID:" $KSW2ID
    mkdir -p $KSW2ID
    echo "Start Date: $d" > $KSW2ID/log.txt
}

run() {
    for dir in $1/*/; do
        # Skip directories that contains "skip" in their name
        is_skip=`echo $dir | grep skip`
        if [[ $is_skip ]]; then
            echo "Skipping $dir" | tee -a $KSW2ID/log.txt
            continue
        fi

        for test in $dir/q*; do
            path=`dirname $test`
            qfile=`basename $test`
            tfile=`echo $qfile | sed 's/q/t/'`
            case=`echo $qfile | sed 's/q//' | sed 's/\.fa//'`
            echo "$path[$case]" | tee -a $KSW2ID/log.txt
            cmd=`echo $KSW2 $path/$qfile $path/$tfile ${@:2}`
            echo "  $cmd" >> $KSW2ID/log.txt
            mkdir -p $KSW2ID/$path
            eval $cmd > $KSW2ID/$path/output$case.txt
        done
    done
}

finalize() {
    d=`date`
    echo "End Date: $d" >> $KSW2ID/log.txt
}

if [ -z  "$1" ]; then
    echo "Usage: $0 <dir> [args]"
    echo "Run a KSW2 Test suite"
    echo ""
    echo "  <dir>   KSW2 Test suite directory"
    echo "  [args]  Arguments to pass to KSW2"
    echo ""
    echo "A test suite is a collection of directories with pairs of (q.fa and t.fa) files"
    echo ""
    echo "     test_suite/small/q0.fa, q1.fa, ... qn.fa, t0.fa, t1.fa, ... tn.fa"
    echo "     test_suite/medium/*"
    echo "     test_suite/large/*"
    echo ""
    echo "This script calls KSW2 for each pair of files in the test suite directory"
    echo "Set the environment variable 'KSW2' so that the script can find ksw2-test."
    echo ""
    echo "To skip a directory, mark it by naming it with 'skip', e.g.,"
    echo ""
    echo "     test_suite/skip_large"
    echo ""
    echo "The output of each run is saved in a directory that is generated using a"
    echo "6 character long ID (a hash of Today's date). Set 'KSW2ID' to override."
    echo ""
    echo "The placement of output files mirrors the test_suite input structure, e.g.,"
    echo ""
    echo "129fc3/log.txt : A summary of all cases run"
    echo "129fc3/small/output0.txt : the output of running ksw2-test using (q0.fa, t0.fa)"
    echo ""
    echo "Example usage:"
    echo ""
    echo "     ./run_test_suite.sh test_suite -t extd2_sse"
    echo ""
    echo "All arguments after the first argument (test_suite) are passed to ksw2_test."
    exit;
fi

init
run $@
finalize
