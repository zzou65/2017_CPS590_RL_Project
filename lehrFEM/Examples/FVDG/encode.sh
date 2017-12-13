#!/bin/sh

# Compress the frames of the .avi movie using MPEG4 or MPEG4 V2 codecs from
# the mplayer.
#
#   Usage: encode source.avi dest.avi

# Copyright 2006-2006 Patrick Meury
# SAM - Seminar for Applied Mathematics
# ETH-Zentrum
# CH-8092 Zurich, Switzerland

#
# Set optimal bitrate
#

bitrate=216000

#
# Set MPEG4 codec options (Only works if divx codec is installed)
#

opt="vbitrate=$bitrate:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:scplx_mask=0.1:tcplx_mask=0.1:naq:trell"
codec="mpeg4"

#
# Set Microsoft MPEG4 V2 codec options
#

# opt="vbitrate=$bitrate:mbd=2:keyint=132:vqblur=1.0:cmp=2:subcmp=2:dia=2:mv0:last_pred=3"
# codec="msmpeg4v2"

#
# Clean temporary files that can interfere with the compression phase
#

rm -f divx2pass.log frameno.avi

#
# Compress the avi movie
#

mencoder -ovc lavc -lavcopts vcodec=$codec:vpass=1:$opt -nosound -o /dev/null $1
mencoder -ovc lavc -lavcopts vcodec=$codec:vpass=2:$opt -nosound -o $2 $1

#
# Cleanup
#

rm -f divx2pass.log
