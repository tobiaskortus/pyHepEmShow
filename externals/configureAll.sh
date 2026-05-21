#!/bin/bash

# Always run relative to the externals/ directory
cd "$(dirname "$0")"

./configure.sh msc_no_ad_fix_no   "-DENABLE_MSC=0 -DENABLE_FLUCTUATION=0 -DENABLE_DEBUG_OUTPUT=1 -DENABLE_AD_FIX=0"
./configure.sh msc_yes_ad_fix_no  "-DENABLE_MSC=1 -DENABLE_FLUCTUATION=1 -DENABLE_DEBUG_OUTPUT=1 -DENABLE_AD_FIX=0"
./configure.sh msc_no_ad_fix_yes  "-DENABLE_MSC=0 -DENABLE_FLUCTUATION=0 -DENABLE_DEBUG_OUTPUT=1 -DENABLE_AD_FIX=1"
./configure.sh msc_yes_ad_fix_yes "-DENABLE_MSC=1 -DENABLE_FLUCTUATION=1 -DENABLE_DEBUG_OUTPUT=1 -DENABLE_AD_FIX=1"
