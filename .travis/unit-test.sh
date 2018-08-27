#!/bin/bash

set -e
errors=0

# Run unit tests
python vcfdistil/vcfdistil_test.py || {
    echo "'python python/vcfdistil/vcfdistil_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E vcfdistil/*.py || {
    echo 'pylint -E vcfdistil/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"
