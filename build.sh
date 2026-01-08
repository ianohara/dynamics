#!/bin/bash
set -e

BUILD_DIR=".build"
RUN_TESTS=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --test|-t)
            RUN_TESTS=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--test|-t]"
            exit 1
            ;;
    esac
done

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake .. 2>&1 | grep -v "^CMake Deprecation Warning" | grep -v "^  Compatibility with CMake" | grep -v "^  Update the VERSION" | grep -v "^CMake Warning (dev)" | grep -v "Policy CMP0148" | grep -v "cmake --help-policy" | grep -v "This warning is for project" | grep -v "^$" || true

make

if [[ $RUN_TESTS -eq 1 ]]; then
    ./tests/tests
fi
