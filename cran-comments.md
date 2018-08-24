## Notes

This is an update of an archived package on CRAN.
The package was archived after CRAN sent messages to the 
previous package maintainer warning about check problems.
The previous maintainer had moved from the university
so his email address did not work any longer.

## Test environments
* local macOS (10.13.6) install, R release version
* ubuntu 14.04 (on Travis CI), R release and devel version
* windows Server 2012 R2 x64 (on AppVeyor), R release version
* win-builder, R release and devel version

## R CMD check results
For macOS, Travis CI, and AppVeyor there were no ERRORs, WARNINGs or NOTEs. 

For win-builder release and devel versions there was one note:

-------------------- 
Maintainer: 'Andrew C. Hooker <andrew.hooker@farmbio.uu.se>'

New submission

Package was archived on CRAN

Possibly mis-spelled words in DESCRIPTION:
  NCA (2:8, 10:49, 12:69)
  PKPD (12:50)
  pharmacodynamic (12:33)
  pharmacokinetic (12:5)

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2018-05-03 as check problems were not
    corrected despite reminders.
----------------------

The first part of the note is just stating that I am the maintainer.  

The second and the last part of the note says that this package was
previously archived on CRAN and that this is a new submission.  
The package was archived after CRAN sent messages to the 
previous package maintainer warning about check problems.
The previous maintainer had moved from the university
so his email address did not work any lnger. 
The check problems have been fixed in this version.  

The third part of the note does not recognize a number of words.
- "NCA": is defined in the description and is a well known acronym, 
         see for example: https://www.ncbi.nlm.nih.gov/pubmed/23007438
- "PKPD": is defined in the description and is a well known acronym, 
         see for example: https://en.wikipedia.org/wiki/PK/PD_models
- "pharmacodynamic" and "pharmacokinetic": are well known terms, 
         see for example: https://en.wikipedia.org/wiki/PK/PD_models

## Downstream dependencies
There are currently no downstream dependencies for this package.
