staging area for the Luminosity Convolution Model codebase - 
TODO[ ] transition to the WIT servers. 

notes for initial pulls, 

either (A) install ROOT 5.28 or (B) alter the CMakeList to reflect the version you WANT to use. Otherwise, the FindROOT.CMake will not be able to phone home. 

install cmake and amend appropraite CMake files inorder to reflect your local environment (as-is it reflects rob/rjrivera/code-dojo/ etc )

PLEASE be sure to branch off the mmaster before you start making changes. 

side-note - TMinuit.h (sp?) doesn't appear to be in 6.08 [version I'm using]


