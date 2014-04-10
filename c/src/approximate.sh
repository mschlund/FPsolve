#for filename in lossy_benchmarks/student-cfg/alphabet-ab/*; do ./newton --lossy --file ${filename} >> "ab-trees".log; done

#for filename in lossy_benchmarks/student-cfg/alphabet-abc/*; do ./newton --lossy --file ${filename} >> "abc-trees".log; done

#echo "alphabet 01 trees..."
#for filename in lossy_benchmarks/student-cfg/alphabet-01/*; do ./newton --lossy --file ${filename} >> "01-trees".log; done

#echo "alphabet ab courcelle..."
#for filename in lossy_benchmarks/student-cfg/alphabet-ab/*; do ./newton --lossyC --file ${filename} >> "ab-courcelle".log; done

#echo "alphabet abc courcelle..."
#for filename in lossy_benchmarks/student-cfg/alphabet-abc/*; do ./newton --lossyC --file ${filename} >> "abc-courcelle".log; done

#echo "alphabet 01 courcelle..."
#for filename in lossy_benchmarks/student-cfg/alphabet-01/*; do ./newton --lossyC --file ${filename} >> "01-courcelle".log; done

#echo "comparing approximations 01..."
#for filename in lossy_benchmarks/student-cfg/alphabet-01/*; do ./newton --lossyComp --file ${filename} --file2 ${filename}; done

#echo "comparing approximations abc..."
#for filename in lossy_benchmarks/student-cfg/alphabet-abc/*; do ./newton --lossyComp --file ${filename} --file2 ${filename}; done

#echo "comparing approximations ab..."
#for filename in lossy_benchmarks/student-cfg/alphabet-ab/*; do ./newton --lossyComp --file ${filename} --file2 ${filename}; done

#echo "comparing languages 01..."
#for filename in lossy_benchmarks/student-cfg/alphabet-01/*; do 
#    for filename2 in lossy_benchmarks/student-cfg/alphabet-01/*; do 
#        echo "files: " ${filename} ${filename2} 
#        ./newton --lossyC --file ${filename} --file2 ${filename2} >> "01-courcelle-comparison.log";
#    done
#done

echo "comparing languages abc..."
for filename in lossy_benchmarks/student-cfg/alphabet-abc/*; do 
    for filename2 in lossy_benchmarks/student-cfg/alphabet-abc/*; do 
        echo "files: " ${filename} ${filename2} 
        ./newton --lossyC --file ${filename} --file2 ${filename2} >> "01-courcelle-comparison.log";
    done
done
