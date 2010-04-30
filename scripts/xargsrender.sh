#Run this script like so ls *.pov | xargs -L1 -P4 xargsrender.sh
if [ ! -e ${@%.pov}.png ]; then 
    povray +H768 +W1024 -D $@
fi