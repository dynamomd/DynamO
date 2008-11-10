for dim in 0 1 2; do
    bzcat output.xml.bz2 | xmlstarlet sel -t -v "//Rij[@dimension='$dim']/Histogram" | gawk 'function acos(x) { return atan2((1.-x^2)^0.5,x) } {if (NF==2) print acos($1), $2}' > rij"$dim"angle.dat

    bzcat output.xml.bz2 | xmlstarlet sel -t -v "//Vij[@dimension='$dim']/Histogram" | gawk 'function acos(x) { return atan2((1.-x^2)^0.5,x) } {if (NF==2) print acos($1), $2}' > vij"$dim"angle.dat
done

echo "set yrange [0:1]; set xrange [-1:1]; set grid polar; set size 1, 0.79; set polar; plot 'rij0angle.dat', 'rij1angle.dat', 'rij2angle.dat', 0.5; pause mouse" | gnuplot

echo "set yrange [0:2]; set xrange [-2:2]; set grid polar; set size 1, 0.79; set polar; plot 'vij0angle.dat', 'vij1angle.dat', 'vij2angle.dat', 0.5; pause mouse" | gnuplot
