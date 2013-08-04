N="10 20 30 40 50 60 70 80 90 100"
eps="0.1 0.2 0.3 0.4 0.5"
#eps="0.02 0.04 0.06 0.08 0.1"

for n in $N
do
  for e in $eps
  do
    python gen_eqs.py $n $e > "c/test/grammars/float-random/float2_${n}_${e/./}.cfg"
  done
done
