#!/bin/bash

ePCC=ePCC

#yes | cp "/mnt/hgfs/D/20210616 PCC/PCC.py" .
pwd=`pwd`
dos2unix -q *
done=done
mkdir ${done}

main=PCC.py
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

echo [${ePCC}] Initial prework...
python ${main} 0 | tee -a ${python_log}

for (( i=0; i<=100; i++ ))
do
echo [${ePCC}] ${i} iteration...

python ${main} ${i} | tee -a ${python_log}
cov=`grep -c -s 'PCC converged!' ${python_log}`

if [ ${cov} -ge 1 ];then
cp ${pwd}/${i}/${i}.chg ${pwd}/conv.chg
break
fi

cd ${i}
dos2unix -q *

echo [${ePCC}] Iteration ${i}, Gaussian calculation...
gaucov=`grep -c -s 'Normal termination' ${i}.log`
if [[ "${gaucov}" -eq 0 ]]; then
g16 < ${i}.gjf 1> ${i}.log 2>&1
else
echo [${ePCC}] Iteration ${i}, Gaussian calculation already completed! 
fi
gaucov=`grep -c -s 'Normal termination' ${i}.log`
if [[ "${gaucov}" -eq 0 ]]; then
echo [ERROR] Iteration ${i}, Gaussian not converged, ${ePCC} exit! 
exit; fi

chgcov=`grep -c -s 'If outputting atom coordinates' ${i}_resp.log`
if [[ "${chgcov}" -eq 0 ]]; then
Multiwfn ${i}.chk < ${resp_input} 1> ${i}_resp.log 2>&1
else
echo [${ePCC}] Iteration ${i}, RESP already fitted! 
fi

chgcov=`grep -c -s 'If outputting atom coordinates' ${i}_resp.log`
if [[ "${chgcov}" -eq 0 ]]; then
echo [ERROR] Iteration ${i}, RESP fitting failure, check log please! ${ePCC} exit! 
exit; fi

cd ${pwd}
done

cd ${pwd}
if [[ "${cov}" -gt 0 ]]; then
python ${main} p
else
echo 'PCC not converged in specified steps! '
exit
fi

for (( j=0; j<=${i}; j++ ))
do
mv ${j} ${done}
done

rm -f *fuse*
rm -rf __pycache__

cp PCC_para.py ${done}/
mv resp.in sym_atom PCC_python.log ${done}
mv *.svg ${done}
cp *CAR ${done}
mv *.chg ${done}
