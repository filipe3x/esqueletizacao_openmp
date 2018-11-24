awk -F, -v RS='\n' '{ for (i=1; i<NF; i++) { print "("i " , " $i")" } print "(" i " , " $NF; print ")"; print " "; }' $1

