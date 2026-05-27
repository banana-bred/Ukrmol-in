#!/bin/bash
#
# Copy results from the current run of the UKRmol+ test suite to the repository.
#
# This is particularly useful when creating the benchmark outputs to be stored
# in the repository. However, given the fact that the numerical results are somewhat
# sensitive to the software stack used, it is useful for regression testing with
# numdiff even during development -- a developer may generate some reference results
# using his/her own software stack and using a known-to-be-correct revision of the
# code and use this script to replace the reference results stored in the repository
# for later use of NUMDIFF.
#
# The script is expected to be executed in the "Testing" directory created by ctest.
#
# For generation of the reference results, only the "serial" subset of the test suite
# should be executed.
#
# Usage:
#
#   > copy-test-suite-reference.sh path/to/ukrmol-in-repo path/to/ukrmol-out-repo
#

if [ -z "$1" ] || [ -z "$2" ]
then
    echo "Usage"
    echo "  > copy-test-suite-reference.sh path/to/ukrmol-in-repo path/to/ukrmol-out-repo"
    exit 1
fi

if [ "$(basename $(pwd))" != "Testing" ]
then
    echo "The current directory is not Testing"
    exit 1
fi

if [ ! -d "$1/tests/suite" ]
then
    echo "The path \"$1/tests/suite\" does not exist"
    exit 1
fi

if [ ! -d "$2/tests/suite" ]
then
    echo "The path \"$2/tests/suite\" does not exist"
    exit 1
fi

for DIR in $(ls -d $1/tests/suite/*/)
do
    tname=$(basename $DIR)
    for par in serial parallel
    do
        tdir=${tname}_${par}
        if [ -d $tdir ]
        then
            echo "UKRmol-in: $tdir"
            cp --verbose \
                $tdir/outputs/target.* \
                $tdir/outputs/scattering.eigs.* \
                $1/tests/suite/$tname/outputs/ \
                2> /dev/null
        fi
    done
done

for DIR in $(ls -d $2/tests/suite/*/)
do
    tname=$(basename $DIR)
    for par in serial parallel
    do
        tdir=${tname}_${par}
        if [ -d $tdir ]
        then
            echo "UKRmol-out: $tdir"
            cp --verbose \
                $tdir/outputs/scattering.eigenp.* \
                $tdir/outputs/scattering.ixsecs.* \
                $tdir/outputs/photoionization.* \
                $tdir/outputs/molecular_data.* \
                $2/tests/suite/$tname/outputs/ \
                2> /dev/null
        fi
    done
done

cp --verbose \
    output.* \
    $1/tests/suite/test_benchmark_outputs/ \
    2> /dev/null
