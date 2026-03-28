#!/usr/bin/perl
# Generated automatically from fxalogdir.pl.in by configure.
$LAPSROOT = $ENV{LAPSROOT};
require "$LAPSROOT/etc/fxa.pm";
$LAPS_LOG_PATH = &Set_logdir'fxa; #'
print "$LAPS_LOG_PATH";
