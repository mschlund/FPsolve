# usage: extract and convert the cfg-analyzer benchmarks using the python script in the same directory (convert_cfg.py)
# put them in the respecive test directory or change the locations below

for f in `find c/test/grammars/cfg-analyzer-benchmarks/benchmarks -name '*.cfg'`
do
  echo $f
  echo "sl-nopt:"
  timeout 15s ./newton -f $f --slset --scc > ${f}_log_slset_nonopt
  if [ $? -ne 0 ]; then
    #!= 0 signals a timeout or that sth else is wrong
    echo "TIMEOUT"
    echo "TIMEOUT: 999111999 ms" >> ${f}_log_slset_nonopt
  fi

  echo "sl-opt:"
  timeout 15s ./newton -f $f --slset --lin-simpl --scc > ${f}_log_slset_simpl
  if [ $? -ne 0 ]; then
    #!= 0 signals a timeout or that sth else is wrong
    echo "TIMEOUT"
    echo "TIMEOUT: 999111999 ms" >> ${f}_log_slset_simpl
  fi

  echo "ml-nopt:"
  timeout 15s ./newton -f $f --mlset --scc > ${f}_log_mlset_nonopt
  if [ $? -ne 0 ]; then
    #!= 0 signals a timeout or that sth else is wrong
    echo "TIMEOUT"
    echo "TIMEOUT: 999111999 ms" >> ${f}_log_mlset_nonopt
  fi

  echo "ml-opt:"
  timeout 15s ./newton -f $f --mlset --vec-simpl --scc > ${f}_log_mlset_simpl
  if [ $? -ne 0 ]; then
    #!= 0 signals a timeout or that sth else is wrong
    echo "TIMEOUT"
    echo "TIMEOUT: 999111999 ms" >> ${f}_log_mlset_simpl
  fi

done

rm runtimes_*_*_*.log

for f in `find c/test/grammars/cfg-analyzer-benchmarks -name *_log_slset_nonopt*`
do
  runtime_ms=`grep -o "[[:digit:]]\{1,\}[[:space:]]ms" $f`
  x=${f##*alphabet-}
  y=${x%%.*}
  z=${y/\//_};
  cfg_file=$z
  echo $cfg_file,${runtime_ms%% ms} >> runtimes_slset_nonopt.log
done

for f in `find c/test/grammars/cfg-analyzer-benchmarks -name *_log_slset_simpl*`
do
  runtime_ms=`grep -o "[[:digit:]]\{1,\}[[:space:]]ms" $f`
  x=${f##*alphabet-}
  y=${x%%.*}
  z=${y/\//_};
  cfg_file=$z
  echo $cfg_file,${runtime_ms%% ms} >> runtimes_slset_simpl.log
done

for f in `find c/test/grammars/cfg-analyzer-benchmarks -name *_log_mlset_simpl*`
do
  runtime_ms=`grep -o "[[:digit:]]\{1,\}[[:space:]]ms" $f`
  x=${f##*alphabet-}
  y=${x%%.*}
  z=${y/\//_};
  cfg_file=$z
  echo $cfg_file,${runtime_ms%% ms} >> runtimes_mlset_simpl.log
done

for f in `find c/test/grammars/cfg-analyzer-benchmarks -name *_log_mlset_nonopt*`
do
  runtime_ms=`grep -o "[[:digit:]]\{1,\}[[:space:]]ms" $f`
  x=${f##*alphabet-}
  y=${x%%.*}
  z=${y/\//_};
  cfg_file=$z
  echo $cfg_file,${runtime_ms%% ms} >> runtimes_mlset_nonopt.log
done



