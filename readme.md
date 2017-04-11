staging area for the Luminosity Convolution Model codebase - 
TODO[ ] transition to the WIT servers. 

notes for initial pulls, 

either (A) install ROOT 5.28 or (B) alter the CMakeList to reflect the version you WANT to use. Otherwise, the FindROOT.CMake will not be able to phone home. 

install cmake and amend appropraite CMake files inorder to reflect your local environment (as-is it reflects rob/rjrivera/code-dojo/ etc )

PLEASE be sure to branch off the mmaster before you start making changes. 

side-note - TMinuit.h (sp?) doesn't appear to be in 6.08 [version I'm using]

TMinuit.h is there, but the old cmake settings aren't finding it properly.  You should be all set now.  Run

```
cmake build/
make
```

That should do it

# Configuration Management Action Items

[ ] pull most recent modifications to CalcH.cc
[ ] provide an image to publicly accessible domain - containing installations of all software prereqs and instructions for branching from master
[ ] investigate the validity and performance/accuracy deltas migrating codebase to GSL.
[ ] refactor codebase for modularity *heavy-duty action-item*


ADDENDUM:
If you get an error regarding failure to load "libCore.so", then your environment variable is not pointing to the newly installed ROOT on your system.

to do so in ubuntu - 

[ ] sudo vim /etc/ld.so.conf.d/<pickaname>.conf
    in the file, write the direct path to the lib "/home/linux/lib"

if you don't know where it is, navigate to home and do a search for libCore.so => find name "libCore.so"

cheers!


