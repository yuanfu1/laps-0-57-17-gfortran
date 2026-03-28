#!/bin/sh

# Run sched.pl primarily for use by parent XML script

SCHEDARGS="$1"  

LAPSINSTALLROOT=$2

LAPS_DATA_ROOT=$3

echo "/usr/bin/perl $LAPSINSTALLROOT/sched.pl $SCHEDARGS $LAPSINSTALLROOT $LAPS_DATA_ROOT"
      /usr/bin/perl $LAPSINSTALLROOT/sched.pl $SCHEDARGS $LAPSINSTALLROOT $LAPS_DATA_ROOT 

touch $LAPS_DATA_ROOT/log/sched.done

