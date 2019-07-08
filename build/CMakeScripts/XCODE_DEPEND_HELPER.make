# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.SQPHOTSTART.Debug:
/Users/xinyiluo/Dropbox/Research/SQPhotstart-LANL/bin/Release/Debug/SQPHOTSTART:
	/bin/rm -f /Users/xinyiluo/Dropbox/Research/SQPhotstart-LANL/bin/Release/Debug/SQPHOTSTART


PostBuild.SQPHOTSTART.Release:
/Users/xinyiluo/Dropbox/Research/SQPhotstart-LANL/bin/Release/Release/SQPHOTSTART:
	/bin/rm -f /Users/xinyiluo/Dropbox/Research/SQPhotstart-LANL/bin/Release/Release/SQPHOTSTART


PostBuild.SQPHOTSTART.MinSizeRel:
/Users/xinyiluo/Dropbox/Research/SQPhotstart-LANL/bin/Release/MinSizeRel/SQPHOTSTART:
	/bin/rm -f /Users/xinyiluo/Dropbox/Research/SQPhotstart-LANL/bin/Release/MinSizeRel/SQPHOTSTART


PostBuild.SQPHOTSTART.RelWithDebInfo:
/Users/xinyiluo/Dropbox/Research/SQPhotstart-LANL/bin/Release/RelWithDebInfo/SQPHOTSTART:
	/bin/rm -f /Users/xinyiluo/Dropbox/Research/SQPhotstart-LANL/bin/Release/RelWithDebInfo/SQPHOTSTART




# For each target create a dummy ruleso the target does not have to exist
