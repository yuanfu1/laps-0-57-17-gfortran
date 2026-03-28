# LAPS Build System Changes for macOS ARM64 with gfortran

This document summarizes all changes made to successfully build LAPS on macOS ARM64 using gfortran and Homebrew libraries.

## Date: November 23, 2025

---

## 1. Configuration Template Changes (PERSISTENT across configure runs)

### File: `src/include/makefile.inc.in`

**Purpose**: Template file used by configure script to generate `src/include/makefile.inc`

#### Changes Made:

1. **C Compiler Flags** (Line 8):
   - **Added**: `-Wno-error=int-to-pointer-cast -Wno-error=pointer-to-int-cast -Wno-error=implicit-function-declaration`
   - **Reason**: Allow legacy C code to compile with modern Clang/GCC strict warnings
   ```makefile
   CFLAGS=@CFLAGS@  @C_OPT@ -Wno-error=int-to-pointer-cast -Wno-error=pointer-to-int-cast -Wno-error=implicit-function-declaration
   ```

2. **Fortran Compiler Flags** (Line 19):
   - **Added**: `-fallow-argument-mismatch`
   - **Reason**: gfortran 10+ requires this flag to allow type/rank mismatches that were permitted by Intel/PGI compilers
   ```makefile
   FFLAGS = @FFLAGS@ -fallow-argument-mismatch $(INC) @F_OPT@
   ```

3. **NetCDF Fortran Library** (Line 63):
   - **Changed**: `-lnetcdf` → `-lnetcdff -lnetcdf`
   - **Reason**: Homebrew separates NetCDF C and Fortran libraries; both needed for Fortran programs
   ```makefile
   OTHERLIBS = $(LIBPATHFLAG)$(NETCDF)/lib -lnetcdff -lnetcdf  @OTHERLIBS@
   ```

---

## 2. Subdirectory Makefile Changes (Static source files - already persistent)

### Files Modified to Clean `*.mod` Files:

Fortran 90/95 generates `.mod` (module interface) files that need cleanup.

1. **`src/mesowave/recurs_iter/Makefile`** (Line 73):
   ```makefile
   clean:
       rm -f $(EXE)  *.o *.mod *~ *#
   ```

2. **`src/humid/Makefile`** (Line 81):
   ```makefile
   clean:
       rm -f $(EXE) $(LIB) *.o *.mod *~ *# *.i
   ```

3. **`src/lib/cloud/Makefile`** (Line 60):
   ```makefile
   clean:
       rm -f $(LIB) $(DEBUGLIB) *.o *.mod *~ *#
   ```

**Note**: Other Makefiles already had `*.mod` in clean targets.

---

## 3. Source Code Fixes (PERSISTENT - part of source tree)

### Fortran Files:

1. **`src/lib/laps_io.f`**
   - Added `implicit none` to 5 wrapper functions
   - Added explicit variable declarations
   - Used temporary 1-element arrays to convert scalars to rank-1 arrays for NetCDF calls
   - Fixed 72-character line limit violations in fixed-form Fortran
   - Fixed variable name typos (i4diff → i4_diff)

2. **`src/lib/read_namelist.f`**
   - Added `implicit none` to first subroutine
   - Carefully ordered declarations before include files containing namelists
   - Added missing local variable declarations

### C Files:

1. **`src/lib/getfilenames_c.c`**
   - Added `#include <string.h>` and `#include <stdlib.h>`
   - Added forward declaration for `nstrncpy()`

2. **`src/lib/regex.c`**
   - Defined `STDC_HEADERS` to ensure stdlib.h is properly included

3. **`src/var/getiofile.c`**
   - Added `#include <unistd.h>` for `read()` and `close()` functions

---

## 4. How to Use These Changes

### For Future Configure Runs:

1. The changes to `src/include/makefile.inc.in` are permanent
2. Simply run configure as normal:
   ```bash
   ./configure --prefix=/path/to/install --with-netcdf=/opt/homebrew
   ```
3. The generated `src/include/makefile.inc` will automatically include the compiler flags

### For Building:

```bash
make clean    # Cleans all .o, .mod, and executables
make          # Builds all libraries and executables
make install  # Installs to configured prefix
```

### For Cleaning:

```bash
make clean      # Remove build artifacts (*.o, *.mod, *.exe)
make realclean  # Clean + clean libraries
make distclean  # Complete cleanup including configure-generated files
```

---

## 5. Key Technical Details

### Why These Changes Were Needed:

1. **gfortran Strictness**:
   - gfortran enforces Fortran standard more strictly than Intel ifort or PGI pgfortran
   - Type and rank mismatches that were warnings in legacy compilers are errors in gfortran 10+
   - Solution: `-fallow-argument-mismatch` flag

2. **Fixed-Form Fortran Constraints**:
   - Code must be in columns 7-72 (66 characters maximum)
   - Long lines were being truncated, cutting variable names mid-token
   - Solution: Carefully split long lines with continuation characters

3. **Implicit Typing Issues**:
   - Fortran's implicit typing (i-n = INTEGER, others = REAL) caused array declarations to be interpreted as function calls
   - Solution: Add `implicit none` and explicit declarations

4. **NetCDF Fortran Library**:
   - Homebrew separates libnetcdf (C) and libnetcdff (Fortran)
   - Both libraries must be linked for Fortran programs using NetCDF
   - Solution: Added `-lnetcdff` before `-lnetcdf`

5. **C Compiler Warnings as Errors**:
   - Modern Clang treats implicit function declarations as errors
   - Legacy C code missing header includes
   - Solution: Add missing headers and -Wno-error flags

---

## 6. Build Success Summary

**26 executables successfully built**, including:
- STMAS3D.exe (4D variational analysis)
- temp.exe, humid.exe (temperature/humidity analysis)
- Various verification and ingest tools
- Data assimilation and grid processing utilities

**All libraries compiled successfully**:
- liblaps.a (78KB)
- libmodules.a, libwind.a, libcloud.a, etc.

---

## 7. Platform Information

- **OS**: macOS ARM64 (Apple Silicon)
- **Compiler**: gfortran (GNU Fortran from Homebrew)
- **C Compiler**: gcc (actually Clang on macOS)
- **Libraries**: Homebrew packages in /opt/homebrew
  - netcdf
  - netcdf-fortran
  - jasper
  - libpng

---

## 8. References

- Fortran Standard: Fixed-form rules (columns 7-72)
- gfortran documentation: Type/rank mismatch handling
- NetCDF Fortran API: Requires array arguments, not scalars
- Homebrew library locations: /opt/homebrew/lib, /opt/homebrew/include

---

*Document prepared: November 23, 2025*
*LAPS Version: 0.57.17*
