#!/bin/bash
../bin/dynamod -m 16 --i2 2 --i1 1 -d 0.14 -C 8 --text
boxlength=$(bzcat config.out.xml.bz2 | xml sel -t -v '//Units/@BoxLength')


echo "Lennard Jones config " > config.gro
bzcat config.out.xml.bz2 | xml sel -t -v '//ParticleData/@N' >> config.gro
bzcat config.out.xml.bz2 \
    | xml sel -t -m '//Pt' -v '@ID' -o ' ' -v 'P/@x' -o ' ' -v 'P/@y' -o ' ' \
    -v 'P/@z' -o ' ' -v 'V/@x' -o ' ' -v 'V/@y' -o ' ' -v 'V/@z' -n \
    | gawk '{if (NF) print ($1+1)"Ar", "Ar1", $1+1, $2+('$boxlength'/2.0),$3+('$boxlength'/2.0),$4+('$boxlength'/2.0),$5+('$boxlength'/2.0),$6+('$boxlength'/2.0),$7+('$boxlength'/2.0)}' \
  >> config.gro
echo "$boxlength $boxlength $boxlength" | gawk '{print $1+0,$2+0,$3+0}' >> config.gro
