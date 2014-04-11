#for filename in lossy_benchmarks/student-cfg/alphabet-ab/*; do ./newton --lossy --file ${filename} >> "ab-trees".log; done

#for filename in lossy_benchmarks/student-cfg/alphabet-abc/*; do ./newton --lossy --file ${filename} >> "abc-trees".log; done

#echo "alphabet 01 trees..."
#for filename in lossy_benchmarks/student-cfg/alphabet-01/*; do ./newton --lossy --file ${filename} >> "01-trees".log; done

#echo "alphabet ab courcelle..."
#for filename in lossy_benchmarks/student-cfg/alphabet-ab/*; do ./newton --lossyC --file ${filename} >> "ab-courcelle".log; done

#echo "alphabet abc courcelle..."
#for filename in lossy_benchmarks/student-cfg/alphabet-abc/*; do ./newton --lossyC --file ${filename} >> "abc-courcelle".log; done
#echo "abc Courcelle..........."
#for filename in lossy_benchmarks/student-cfg/alphabet-abc/*; do ./newton --lossyC --file ${filename}; done
#echo -e "\n\n\nab Courcelle............"
#for filename in lossy_benchmarks/student-cfg/alphabet-ab/*; do ./newton --lossyC --file ${filename}; done
#echo -e "\n\n\n01 Courcelle............"
#for filename in lossy_benchmarks/student-cfg/alphabet-01/*; do ./newton --lossyC --file ${filename}; done
#echo -e "\n\n\na Courcelle............"
#for filename in lossy_benchmarks/student-cfg/alphabet-a/*; do ./newton --lossyC --file ${filename}; done

echo "abc Derivation Trees..........."
for filename in lossy_benchmarks/student-cfg/alphabet-abc/*; do ./newton --lossy --file ${filename}; done
echo -e "\n\n\nab Derivation Trees............"
for filename in lossy_benchmarks/student-cfg/alphabet-ab/*; do ./newton --lossy --file ${filename}; done
echo -e "\n\n\n01 Derivation Trees............"
for filename in lossy_benchmarks/student-cfg/alphabet-01/*; do ./newton --lossy --file ${filename}; done
echo -e "\n\n\na Derivation Trees............"
for filename in lossy_benchmarks/student-cfg/alphabet-a/*; do ./newton --lossyC --file ${filename}; done

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
#    for filename2 in lossy_benchmarks/student-cf g/alphabet-01/*; do 
#        echo "files: " ${filename} ${filename2} 
#        ./newton --lossyC --file ${filename} --file2 ${filename2} >> "01-courcelle-comparison.log";
#    done
#done

#path="lossy_benchmarks/student-cfg/alphabet-abc/student-"
#echo "comparing student languages abc, iterating over (i,j)"
#for i in {1..102}; do
#    for ((j=i+1;j<=103;j++)); do
#        echo "files: " ${path}${i}.g ${path}${j}.g
#        echo -n ${i},${j}": " >> "abc-courcelle-comparison.log"
#        ./newton --lossyC --file ${path}${i}.g --file2 ${path}${j}.g >> "abc-courcelle-comparison.log";
#    done
#done

#path="lossy_benchmarks/student-cfg/alphabet-abc/student-"
#pathA="lossy_benchmarks/answer-cfg/alphabet-abc/answer-"
#echo "comparing student languages with answer languages abc, iterating over (i,j)"
#for i in {1..103}; do
#    for j in {1..3}; do
#        echo "files: " ${path}${i}.g ${pathA}${j}.g
#        echo -n ${i},${j}": " >> "abc-courcelle-student-answer-comparison.log"
#        ./newton --lossyC --file ${path}${i}.g --file2 ${pathA}${j}.g >> "abc-courcelle-student-answer-comparison.log";
#    done
#done

path="lossy_benchmarks/student-cfg/alphabet-01/student-"
pathA="lossy_benchmarks/answer-cfg/alphabet-01/answer-"
#echo "comparing student languages with answer languages 01, iterating over (i,j)"
#for i in {1..350}; do
#    for j in {1..13}; do
#        echo "files: " ${path}${i}.g ${pathA}${j}.g
#        echo -n ${i},${j}": " >> "01-courcelle-student-answer-comparison.log"
#        ./newton --lossyC --file ${path}${i}.g --file2 ${pathA}${j}.g >> "01-courcelle-student-answer-comparison.log";
#    done
#done

#for j in {1..13}; do
#    ./newton --lossyC --file ${pathA}${j}.g;
#done 

#echo "comparing languages abc..."
#for filename in lossy_benchmarks/student-cfg/alphabet-abc/*; do 
#    for filename2 in lossy_benchmarks/student-cfg/alphabet-abc/*; do 
#        echo "files: " ${filename} ${filename2} 
#        ./newton --lossyC --file ${filename} --file2 ${filename2} >> "abc-courcelle-comparison.log";
#    done
#done
