#!/bin/sh
make -C /Users/xinyiluo/Dropbox/Research/SQPhotstart-LANL/build -f /Users/xinyiluo/Dropbox/Research/SQPhotstart-LANL/build/CMakeScripts/ZERO_CHECK_cmakeRulesBuildPhase.make$CONFIGURATION OBJDIR=$(basename "$OBJECT_FILE_DIR_normal") all
