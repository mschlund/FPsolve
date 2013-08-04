BENCH_PREFIX='float_'

for f in `find c/test/grammars/float-random -name ${BENCH_PREFIX}'*.cfg'`
do
  echo $f
  echo "float:"
  timeout 3m ./c/src/newton -f $f --float -i 10 > ${f}_log_float_noscc
  if [ $? -ne 0 ]; then
    #!= 0 signals a timeout or that sth else is wrong
    echo "TIMEOUT"
    echo "TIMEOUT: 999111999 ms" >> ${f}_log_float_noscc
  fi
done

rm runtimes_${BENCH_PREFIX}.log

for f in `find c/test/grammars/float-random -name ${BENCH_PREFIX}'*_log_float_noscc*'`
do
  runtime_ms=`grep -o "[[:digit:]]\{1,\}[[:space:]]ms" $f`
  x=${f##*${BENCH_PREFIX}_}
  y=${x%%.*}
  z=${y/\//_};
  cfg_file=$z
  echo $cfg_file,${runtime_ms%% ms} >> runtimes_${BENCH_PREFIX}_noscc.log
done

sort -g runtimes_${BENCH_PREFIX}_noscc.log > rt.tmp
mv rt.tmp runtimes_${BENCH_PREFIX}_noscc.log
