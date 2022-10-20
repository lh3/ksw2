#!/bin/bash

pass=0
total=0
error=0

summary() {
    for dir in $1//*/*/; do
        for test in ${dir}/output*; do 
            test=`echo $test | sed -E 's|/+|/|g'`
            ref=`echo $test | sed -E "s|$1|$2/|"`
            test_count=`cat $ref | wc -l`
            test_errors=`diff $test $ref -y --suppress-common-lines | wc -l`
            difference=`diff $test $ref`
            test_pass=$((test_count-test_errors))
            if [[ $test_pass != $test_count ]]; then
                printf "$test \u1b[31m F [$test_pass / $test_count] \033[m \n"
                error=$((error+1))
            else
                printf "$test \u1b[32m T [$test_pass / $test_count] \033[m \n"
                pass=$((pass+1))
            fi

            total=$((total+1))
        done
    done
}

details() {
    for dir in $1//*/*/; do
        for test in $dir/output*; do 
            test=`echo $test | sed -E 's|/+|/|g'`
            ref=`echo $test | sed -E "s|$1|$2/|"`
            difference=`diff $test $ref`

            if [[ $difference ]]; then
                printf "Error in: $test:"
                diff --unchanged-line-format="" --old-line-format="" --new-line-format="%dn, " \
                     $test $ref
                printf "\n"
            fi
        done
    done
}

statusbar() {
        if [[ $pass != $total ]]; then
            printf "\u1b[31m================================================="
            printf "\u1b[31m $error failed\033[m"
            printf ", "
            printf "\u1b[32m $pass passed\033[m "
            printf "\u1b[31m================================================="
            printf "\033[m\n"
        else
            printf "\u1b[32m================================================="
            printf "\u1b[32m $pass passed\033[m "
            printf "\u1b[32m================================================="
            printf "\033[m\n"
        fi
}

if [ -z  "$1" ];
then
    echo "Usage: $0 <test> <reference> [opt]"
    echo "Compare KSW2 Test suite against reference"
    echo ""
    echo "  test:         Test suite to compare against reference"
    echo "  reference:    Reference test suite"
    echo "  -d:           Show differences (if any)"
    exit;
fi

summary $1 $2
if [ ${3-d} ==  "-d" ]; then
    details $1 $2
fi
statusbar
