#!/bin/bash

ePCC=ePCC

#yes | cp "/mnt/hgfs/D/20210616 PCC/PCC.py" .
pwd=`pwd`
dos2unix *
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

cd 0
echo [${ePCC}] Initial Gaussian calculation...
g16 < 0.gjf 1> 0.log 2>&1
gaucov=`grep -c "Normal termination" 0.log`
if [ ${gaucov} -eq 0 ]; then exit; fi
echo [${ePCC}] Initial Multiwfn calculation...
Multiwfn 0.chk < ${resp_input}
cd ${pwd}

for (( i=1; i<=100; i++ ))
do
echo [${ePCC}] ${i} iteration...

if [ -f "${i}/${i}.chg" ]; then
break
fi

python ${main} ${i} | tee -a ${python_log}
cov=`grep -c "PCC converged!" ${python_log}`

if [ ${cov} -ge 1 ];then
cp ${pwd}/${i}/${i}.chg ${pwd}/conv.chg
break
fi

cd ${i}
dos2unix *

echo [${ePCC}] Iteration ${i}, Gaussian calculation...
g16 < ${i}.gjf 1> ${i}.log 2>&1
gaucov=`grep -c "Normal termination" ${i}.log`
if [ ${gaucov} -eq 0 ]; then exit; fi

Multiwfn ${i}.chk < ${resp_input}
cd ${pwd}
done

cd ${pwd}
if [ ${cov} -ge 1 ]; then
python ${main} p
else
echo 'PCC not converged! '
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
