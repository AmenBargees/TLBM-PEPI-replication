#!/bin/csh -f

set defaultcodename = none

set homedir = $homedir

if(($#argv == 1 && $1 == "-h") || $#argv > 1) then
    if($#argv > 1) then
        echo "--------------------------------------------------------------------------------"
        echo "ERROR: Incorrect usage of make.s. See documentation below for correct usage"
    endif
    echo "--------------------------------------------------------------------------------"
    echo "Compile GTXLBM.c & create the shared object file required to run GTXLBM_run.py"
    echo "--------------------------------------------------------------------------------"
    echo "Usage: make.s [codename]"
    echo " "
    echo "Note: The optional argument codename is the codename for your chips which"
    echo "      defaults to "$defaultcodename" if unspecified. The correct codename"
    echo "      for your CPU is required to activate the AVX-512 vector units & achieve"
    echo '      the highest performance for the "code = c" input parameter for GTXLBM'
    echo '      To verify the correct codename of your cpu, do "lscpu | grep "Model name"'
    echo "      which will tell you what Intel chip you are on such as:"
    echo " "
    echo "Model name:                           Intel (R) Xeon (R) CPU Max 9480"
    echo " "
    echo '   Then google this chip (ie. google Intel Max 9480) to find the "Code Name" for this chip'
    echo "   and look at https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html after -march=cpu-type to"
    echo "   determine the correct codename (cpu-type) of the cpu for the compiler to activate the AVX-512"
    echo "   vector units (ie. sapphirerapids is the codename of the above chip)"
    echo "--------------------------------------------------------------------------------"
    echo "Examples: make.s"
    echo "          make.s sapphirerapids"
    echo "--------------------------------------------------------------------------------"
    exit
endif

echo "--------------------------------------------------------------------------------"
if($#argv == 0) then
    set codename = $defaultcodename
    echo "The CPU chip codename was unspecified"
    echo "--------------------------------------------------------------------------------"
    echo "Note: specification of the codename will activate AVX-512"
    echo "      vector units on Intel chips and optimize performance"
    echo "--------------------------------------------------------------------------------"
    echo "For instructions on how to obtain the chip codename on your computer type make.s -h"
    echo "--------------------------------------------------------------------------------"
else if($#argv == 1) then
    set codename = $1
    echo "The CPU chip codename was specified as $codename"
    echo "--------------------------------------------------------------------------------"
else
    echo 'Incorrect usage. Useage is make.s or "make.s codename" where codename is the CPU chip code name'
    echo "--------------------------------------------------------------------------------"
    exit
endif

pushd $homedir >& /dev/null

echo -n "Compiling GTXLBM.c ... "
if($codename == none) then
    gcc -c -w -fPIC -Ofast GTXLBM.c -lm
else
    gcc -c -w -fPIC -march=$codename -Ofast GTXLBM.c -lm
endif
echo "GTXLBM.o  created"
echo "--------------------------------------------------------------------------------"
echo -n "Compiling GTXLBM.o ... "
gcc -shared GTXLBM.o -lm -o GTXLBM.so
echo "GTXLBM.so created"
rm GTXLBM.o
echo -n "Removing  GTXLBM.o ..."
echo ".......... removed"
echo "--------------------------------------------------------------------------------"

popd >& /dev/null
