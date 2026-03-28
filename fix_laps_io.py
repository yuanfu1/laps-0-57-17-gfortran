#!/usr/bin/env python3
"""
Fix laps_io.f to use temporary arrays instead of scalars when calling READ_LAPS_DATA
"""

# Read the file
with open('src/lib/laps_io.f', 'r') as f:
    lines = f.readlines()

# Function to insert array copy code before and after a READ_LAPS_DATA call
def insert_array_copy(lines, call_line_idx):
    # Find the end of the CALL statement
    i = call_line_idx
    while i < len(lines) and 'ISTATUS)' not in lines[i]:
        i += 1
    
    if i >= len(lines):
        return lines
    
    # Insert copy-back code after the call
    indent = '            '
    copy_back = [
        '\n',
        indent + '!           Copy arrays back to scalars\n',
        indent + 'var_2d = var_2d_arr(1)\n',
        indent + 'lvl_2d = lvl_2d_arr(1)\n',
        indent + 'lvl_coord_2d = lvl_coord_2d_arr(1)\n',
        indent + 'units_2d = units_2d_arr(1)\n',
        indent + 'comment_2d = comment_2d_arr(1)\n',
        '\n',
    ]
    
    # Insert copy-to code before the call
    copy_to = [
        '\n',
        indent + '!           Copy scalars to arrays for READ_LAPS_DATA\n',
        indent + 'var_2d_arr(1) = var_2d\n',
        indent + 'lvl_2d_arr(1) = lvl_2d\n',
        indent + 'lvl_coord_2d_arr(1) = lvl_coord_2d\n',
        indent + 'units_2d_arr(1) = units_2d\n',
        indent + 'comment_2d_arr(1) = comment_2d\n',
        '\n',
    ]
    
    # Replace scalar arguments with array arguments in the CALL
    for j in range(call_line_idx, i+1):
        lines[j] = lines[j].replace('VAR_2D,', 'VAR_2D_arr,')
        lines[j] = lines[j].replace('LVL_2D,', 'LVL_2D_arr,')
        lines[j] = lines[j].replace('LVL_COORD_2D,', 'LVL_COORD_2D_arr,')
        lines[j] = lines[j].replace('UNITS_2D,', 'UNITS_2D_arr,')
        lines[j] = lines[j].replace('COMMENT_2D,', 'COMMENT_2D_arr,')
        # Handle bracket syntax too
        lines[j] = lines[j].replace('[VAR_2D]', 'VAR_2D_arr')
        lines[j] = lines[j].replace('[LVL_2D]', 'LVL_2D_arr')
        lines[j] = lines[j].replace('[LVL_COORD_2D]', 'LVL_COORD_2D_arr')
        lines[j] = lines[j].replace('[UNITS_2D]', 'UNITS_2D_arr')
        lines[j] = lines[j].replace('[COMMENT_2D]', 'COMMENT_2D_arr')
    
    # Insert the code
    lines[call_line_idx:call_line_idx] = copy_to
    lines[i+1+len(copy_to):i+1+len(copy_to)] = copy_back
    
    return lines

# Find lines with CALL READ_LAPS_DATA in 2D wrapper functions
# These are in get_laps_2d (around line 66), get_laps_2dgrid (around line 290), get_laps_2dvar (around line 383)

call_lines = []
for i, line in enumerate(lines):
    if 'CALL READ_LAPS_DATA' in line:
        # Check if this is in a 2D wrapper function (kdim=1, uses VAR_2D scalars)
        # Look at the next few lines for evidence
        context = ''.join(lines[i:min(i+4, len(lines))])
        # Check for 1,1 pattern (first arg is the field count, second is kdim)
        # and VAR_2D (not VAR_3D)
        if 'VAR_2D' in context and ('1,1,' in context or ', 1, 1,' in context or ',1,1 ' in context or ' 1,1,' in context):
            call_lines.append(i)

print(f"Found {len(call_lines)} calls to fix at lines: {[i+1 for i in call_lines]}")

# Process in reverse order so line numbers don't shift
for call_line in reversed(call_lines):
    lines = insert_array_copy(lines, call_line)
    print(f"Fixed call at line {call_line+1}")

# Write the fixed file
with open('src/lib/laps_io.f', 'w') as f:
    f.writelines(lines)

print("Done!")
