#!/bin/csh -f

if($#argv == 1 && $1 == "-h") then
    echo "-------------------------------------------------------------------"
    echo "Install GTXLBM, compilation of C code & specification of GTXLBM homedir"
    echo "-------------------------------------------------------------------"
    echo "Usage: install.s"
    echo "-------------------------------------------------------------------"
    exit
endif

echo set  homedir = `pwd`

echo "================================================================================"
echo "Installing GTXLBM in $homedir ..."

source $homedir/GTXLBM.cshrc

# Compile GTXLBM.c and make the shared object file GTXLBM.so required by GTXLBM_run.py

make.s

echo 'Appending "source $homedir/GTXLBM.cshrc" to your .cshrc file'

set date = `get_date`

echo ""                                                                             >> $HOME/.cshrc
echo "#===========================================================================" >> $HOME/.cshrc
echo "# The following line was added by installation of GTXLBM on $date"            >> $HOME/.cshrc
echo "#---------------------------------------------------------------------------" >> $HOME/.cshrc
echo "source $homedir/GTXLBM.cshrc"                                                 >> $HOME/.cshrc
echo "#===========================================================================" >> $HOME/.cshrc

source $homedir/GTXLBM.cshrc

echo "--------------------------------------------------------------------------------"
echo "Setting GTXLBM home directory to $homedir"

echo $homedir > .homedir

echo "================================================================================"
echo "GTXLBM installation complete in $homedir, GTXLBM_run.py is ready to run"
echo "--------------------------------------------------------------------------------"
echo "GTXLBM aliases will be available and scripts will be in your PATH when you"
echo 'enter new c-shell which can be done by typing "csh". If you are already in'
echo "a c-shell, you can type "source ~/.cshrc". It is recommended to set your"
echo "shell to csh so you will always have access to GTXLBM aliases and scripts."
echo "================================================================================"
