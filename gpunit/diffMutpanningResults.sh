# Copyright (c) 2003-2017 Broad Institute, Inc., Massachusetts Institute of Technology, and Regents of the University of California.  All rights reserved.
#!/bin/sh
execDir=`dirname $0`

zip1=$1
zip2=$2

base1=`basename $1`
base2=`basename $2`
bare1=${base1%%.zip}
bare2=${base2%%.zip}
diffDir1=`mktemp -d $bare1.XXXXXX`
diffDir2=`mktemp -d $bare2.XXXXXX`

tar -xf $zip1 -C $diffDir1
tar -xf $zip2 -C $diffDir2

# the files won't be identical but the top gene should be
cut -f1 $diffDir1/SignificanceSkin.txt | head -2 > $diffDir1/topgenes.txt
cut -f1 $diffDir2/SignificanceSkin.txt | head -2 > $diffDir2/topgenes.txt

# Diff only selected files out of the ZIP
diff --strip-trailing-cr -q $diffDir1/topgenes.txt $diffDir2/topgenes.txt
status=$?


rm -rf $diffDir1 $diffDir2
exit $status
