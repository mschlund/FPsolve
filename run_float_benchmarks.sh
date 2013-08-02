for f in `find c/test/grammars/float-random -name '*.cfg'`
do
  echo $f
  echo "float:"
  timeout 3m ./c/src/newton -f $f --float --scc > ${f}_log_float
  if [ $? -ne 0 ]; then
    #!= 0 signals a timeout or that sth else is wrong
    echo "TIMEOUT"
    echo "TIMEOUT: 999111999 ms" >> ${f}_log_float
  fi
done

rm runtimes_float.log

for f in `find c/test/grammars/float-random -name *_log_float*`
do
  runtime_ms=`grep -o "[[:digit:]]\{1,\}[[:space:]]ms" $f`
  x=${f##*float_}
  y=${x%%.*}
  z=${y/\//_};
  cfg_file=$z
  echo $cfg_file,${runtime_ms%% ms} >> runtimes_float.log
done

sort -g runtimes_float.log > rt.tmp
mv rt.tmp runtimes_float.log
