# LAPS Configuration Fix for macOS

## Problem
The LAPS configure script fails to properly detect type sizes on macOS, setting `SIZEOF_INT`, `SIZEOF_LONG`, and `SIZEOF_SHORT` to 0 in `src/include/config.h`. This causes compilation to fail because `static_routines.c` cannot define the `fint4` type.

## Solution
Two scripts have been created to automatically fix this issue:

### 1. fix_config.sh
Fixes the SIZEOF definitions in config.h after configure runs:
- Sets SIZEOF_INT = 4
- Sets SIZEOF_LONG = 8
- Sets SIZEOF_SHORT = 2

Usage:
```bash
./fix_config.sh
```

### 2. configure_and_fix.sh (Recommended)
Wrapper script that runs configure and automatically applies the fix:

Usage:
```bash
./configure_and_fix.sh [configure options]
```

This is the recommended way to configure LAPS on macOS. It:
1. Runs configure with FC=gfortran F90=gfortran
2. Automatically applies the SIZEOF fixes
3. Reports success/failure

## Files Modified

### Source Code Fixes
The following Fortran files have been fixed for gfortran compatibility:
- `src/lib/gridgen_utils.f` - BOZ literals wrapped with INT()
- `src/lib/process_goes_snd.f` - endian4 subroutine uses transfer()
- `src/lib/get_poes_data.f` - NetCDF scalar variables converted to arrays
- `src/lib/get_rtamps_data.f` - NetCDF scalar variables converted to arrays
- `src/lib/read_rtamps_data.f` - NetCDF scalar variables converted to arrays
- `src/lib/get_radiometer_data.f` - NetCDF scalar variables converted to arrays
- `src/lib/read_local_cwb.f` - Removed undeclared variable initializations
- `src/lib/lvd_sat_ingest.f` - Fixed namelist declaration order

### Build System Fixes
- `configure.in` - Removed `-fd-lines-as-comments` flag (line 342)
  - **Note:** After modifying configure.in, run `autoconf` to regenerate the configure script

## Quick Start

To configure and build LAPS on macOS:

```bash
# Configure with automatic fix
./configure_and_fix.sh

# Build
make
```

## Why This is Needed

The autoconf `AC_CHECK_SIZEOF` macro in `configure.in` (lines 455-457) fails to detect type sizes on modern macOS systems. This appears to be a known issue with older autoconf scripts. The fix scripts provide a reliable workaround by setting the correct values after configure completes.

## Verification

After running configure_and_fix.sh, you can verify the fix:

```bash
grep SIZEOF_ src/include/config.h
```

Should show:
```c
#define SIZEOF_INT 4
#define SIZEOF_LONG 8
#define SIZEOF_SHORT 2
```
