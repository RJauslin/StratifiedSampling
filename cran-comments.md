## Test environments
* Windows 11, R-devel (local)

## R CMD check results and check_win_devel

- There were no ERRORs, WARNINGs or NOTEs

- devtools::check_mac_release
  * checking installed package size ... NOTE
  installed size is  6.3Mb
  sub-directories of 1Mb or more:
    libs   6.0Mb

## rhub

- Windows Servers 2022
  * checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
  'lastMiKTeXException'
  

- Fedora Linux, R-devel, clang, gfortran
 
 * checking installed package size ... NOTE
    installed size is  5.5Mb
    sub-directories of 1Mb or more:
      libs   5.1Mb

  * checking HTML version of manual ... NOTE
      Skipping checking HTML validation: no command 'tidy' found
      
It seems this warning could be supressed by setting _R_CHECK_RD_VALIDATE_RD2HTML_ to false, but apparently that just turns off HTML validation, which happens anyways.
  
- Ubuntu Linux 20.04.1 LTS, R-release, GCC 
- Debian Linux, R-devel, GCC ASAN/UBSAN
 
 * run a PREPERROR, this seems related to rhub issue.