# set up environment variables
# the ${VAR:+:} part adds a double colon only if VAR is not empty
export PATH="/home/pub/B/altel_acts/INSTALL/bin${PATH:+:}${PATH}"
export LD_LIBRARY_PATH="/home/pub/B/altel_acts/INSTALL/lib64${LD_LIBRARY_PATH:+:}${LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="/home/pub/B/altel_acts/INSTALL/lib64${DYLD_LIBRARY_PATH:+:}${DYLD_LIBRARY_PATH}"
