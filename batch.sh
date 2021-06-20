#!/bin/bash

pwd=`pwd`
dos2unix *
rm -f *fuse*

python PCC.py 0

cd 0
g16 < 0.gjf > 0.log
Multiwfn <<EOF
0.chk
7
18
5
1
sym_atom
1
y
0
0
q
EOF
cd ${pwd}

for (( i=0; i<=10; i++ ))
do
echo ${i}
python PCC.py ${i} | tee ${i}_pcc.log
cd ${i}
dos2unix *
g16 < ${i}.gjf > ${i}.log
Multiwfn <<EOF
0.chk
7
18
5
1
sym_atom
1
y
0
0
q
EOF
cd ${pwd}
done
