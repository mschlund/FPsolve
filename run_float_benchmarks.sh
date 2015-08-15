BENCH_PREFIX='float_'

for f in `find c/test/grammars/float-random -name ${BENCH_PREFIX}'*.cfg'`
do
  for i in {1..5}
  do
    echo $f
    echo "float:"
#    outf=${f}_log_float_noscc_${i}
    outf=${f}_log_float_${i}
    echo $outf
#    echo `timeout 5m ./c/src/newton -f $f --float -i 10 | grep "Multiplication\|Addition\|Star" | wc -l`' ms' > $outf
    timeout 5m ./c/src/fpsolve -f $f --float -i 10 -s newtonNumeric > $outf
    if [ $? -ne 0 ]; then
      #!= 0 signals a timeout or that sth else is wrong
      echo "TIMEOUT"
      echo "TIMEOUT: 999111999 ms" >> ${outf}
    fi
  done
done

#rm runtimes_${BENCH_PREFIX}.log

for i in {1..5}
do
  for f in `find c/test/grammars/float-random -name ${BENCH_PREFIX}*log_float_${i}`
  do
    runtime_ms=`grep -o "[[:digit:]]\{1,\}[[:space:]]ms" $f`
    x=${f#*${BENCH_PREFIX}}
    y=${x%%.*}
    z=${y/\//_};
    cfg_file=$z
    echo $cfg_file,${runtime_ms%% ms} >> runtimes_${BENCH_PREFIX}_${i}.log
  done
  sort -g runtimes_${BENCH_PREFIX}_${i}.log > rt.tmp
  mv rt.tmp runtimes_${BENCH_PREFIX}_${i}.log
done

