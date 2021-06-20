#!/bin/bash

pwd=`pwd`
dos2unix *

python_log=${pwd}/PCC_python.log
rm -f ${python_log}
resp_input=${pwd}/resp.in
rm -f ${resp_input}
rm -f sym_atom

cat > ${resp_input} <<EOF
7
18
5
1
${pwd}/sym_atom
1
y
0
0
q

EOF

python PCC.py 0 | tee -a ${python_log}

cd 0
g16 < 0.gjf | tee 0.log
Multiwfn 0.chk < ${resp_input}
cd ${pwd}

for (( i=1; i<=100; i++ ))
do
echo ${i}

if [ -f "${i}/${i}.chg" ]; then
break
fi

python PCC.py ${i} | tee -a ${python_log}
cov=`grep -c "PCC converged!" ${python_log}`

if [ ${cov} -eq 1 ];then
break
fi

cd ${i}
dos2unix *

#last_i=`expr ${i} - 1`
#cp ${pwd}/${last_i}/${last_i}.chk .
g16 < ${i}.gjf | tee ${i}.log
Multiwfn ${i}.chk < ${resp_input}
cd ${pwd}
done

if [ ${cov} -eq 1 ]; then
python PCC.py p
else
echo 'PCC not converged! '
fi

rm -f *fuse*
rm -rf __pycache__
