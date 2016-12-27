#!/bin/sh

if [ ! -d msr-safe ] ; then
    git clone https://github.com/LLNL/msr-safe.git
fi
cd msr-safe
if [ ! -f msr-safe.ko ] ; then
    make
fi
sudo insmod msr-safe.ko
sudo sh -c 'cat whitelists/wl_062D > /dev/cpu/msr_whitelist'

sudo chmod o+rw /dev/cpu/*/msr_safe

echo done

