#!/bin/csh

# LAPS DATA ROOT is first argument

# Path to MSS directory is second argument

umask 002

setenv LAPS_DATA_ROOT $1
setenv MSSPATH $2

#setenv LAPSINSTALLROOT /usr/nfs/lapb/builds/laps

setenv YYDDDHHMM `/usr/bin/perl /Users/xiey/developments/da/laps/lapsinstallroot/etc/sched_sys.pl -d 1.00 -D $LAPS_DATA_ROOT`
setenv MONTH     `/usr/bin/perl /Users/xiey/developments/da/laps/lapsinstallroot/etc/sched_sys.pl -d 1.00 -f mm -D $LAPS_DATA_ROOT`
setenv YYYYMMDD_HHMM `/usr/bin/perl /Users/xiey/developments/da/laps/lapsinstallroot/etc/sched_sys.pl -d 1.00 -f yyyymmdd_hhmm -D $LAPS_DATA_ROOT`

echo "laps time $YYDDDHHMM"

echo "a13time $YYYYMMDD_HHMM"

echo "month = $MONTH"

setenv YYYY `echo $YYYYMMDD_HHMM  | cut -c1-4`

setenv DATE `echo $YYYYMMDD_HHMM  | cut -c7-8`

echo "MSSFULLPATH = $MSSPATH/$YYYY/$MONTH"

cd $LAPS_DATA_ROOT

echo "/Users/xiey/developments/da/laps/lapsinstallroot/etc/tarlapstime.sh $LAPS_DATA_ROOT $YYDDDHHMM $YYYY $MONTH$DATE lga noexpand nostatic 48 112"
      /Users/xiey/developments/da/laps/lapsinstallroot/etc/tarlapstime.sh $LAPS_DATA_ROOT $YYDDDHHMM $YYYY $MONTH$DATE lga noexpand nostatic 48 112

mssMkdir -p $MSSPATH/$YYYY/$MONTH

mssPut laps_$YYDDDHHMM.tar* $MSSPATH/$YYYY/$MONTH

rm -f laps_$YYDDDHHMM.tar*

