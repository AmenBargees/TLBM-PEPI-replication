#!/bin/csh -f

if(($#argv == 1 && $1 == "-h") || $#argv > 2) then
    echo "-------------------------------------------------------------------"
    echo "Get your computer's IP address"
    echo "-------------------------------------------------------------------"
    exit
endif

set ip_address = `ifconfig | sed -En 's/127.0.0.1//;s/.*inet (addr:)?(([0-9]*\.){3}[0-9]*).*/\2/p'`

# Alternative methods to get the IP address as belo
# ip addr
# ifconfig ifconfig.me

echo "-------------------------------------------------------------------"
echo "IP address: $ip_address"
echo "-------------------------------------------------------------------"
echo "If your IP address is blank, you are not connected to the internet"
echo "-------------------------------------------------------------------"
