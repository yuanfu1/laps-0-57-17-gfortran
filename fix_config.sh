#!/bin/bash
# Fix SIZEOF definitions in config.h after configure runs
# This is needed because AC_CHECK_SIZEOF fails on macOS

CONFIG_H="src/include/config.h"

if [ ! -f "$CONFIG_H" ]; then
    echo "Error: $CONFIG_H not found"
    exit 1
fi

# Use sed to replace the SIZEOF definitions
sed -i.bak \
    -e 's/^#define SIZEOF_INT 0$/#define SIZEOF_INT 4/' \
    -e 's/^#define SIZEOF_LONG 0$/#define SIZEOF_LONG 8/' \
    -e 's/^#define SIZEOF_SHORT 0$/#define SIZEOF_SHORT 2/' \
    "$CONFIG_H"

echo "Fixed SIZEOF definitions in $CONFIG_H"
echo "  SIZEOF_INT = 4"
echo "  SIZEOF_LONG = 8"  
echo "  SIZEOF_SHORT = 2"
