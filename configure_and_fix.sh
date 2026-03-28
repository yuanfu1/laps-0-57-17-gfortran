#!/bin/bash
# Wrapper script to run configure and fix SIZEOF definitions
# Usage: ./configure_and_fix.sh [configure options]

echo "Running configure..."
FC=gfortran F90=gfortran ./configure "$@"
CONFIGURE_STATUS=$?

if [ $CONFIGURE_STATUS -ne 0 ]; then
    echo "Configure failed with status $CONFIGURE_STATUS"
    exit $CONFIGURE_STATUS
fi

echo ""
echo "Applying SIZEOF fixes..."
./fix_config.sh

echo ""
echo "Configuration complete!"
