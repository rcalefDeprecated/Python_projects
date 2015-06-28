#/bin/sh
./count-kmers -o0 < mark_f14_1.seqs > counts_0_f_14_1
./count-kmers -o1 < mark_f14_1.seqs > counts_o1_f_14_1
./count-kmers -o0 < mark_f14_2.seqs > counts_0_f_14_2 
./count-kmers -o1 < mark_f14_2.seqs > counts_o1_f_14_2
./count-kmers -o0 --alphabet=ACDEFGHIKLMNPQRSTVWY < mark_f14_1.seqs > counts_0_f_14_1_constrained
./count-kmers -o1 --alphabet=ACDEFGHIKLMNPQRSTVWY < mark_f14_1.seqs > counts_o1_f_14_1_constrained
./count-kmers -o0 --alphabet=ACDEFGHIKLMNPQRSTVWY < mark_f14_2.seqs > counts_0_f_14_2_constrained
./count-kmers -o1 --alphabet=ACDEFGHIKLMNPQRSTVWY < mark_f14_2.seqs > counts_o1_f_14_2_constrained

echo Train mark_f14_1.seqs Order 0 > output.test
./encoding-cost mark_f14_2.seqs  mark_f14_1.seqs < counts_0_f_14_1 >> output.test
echo >> output.test
echo Train mark_f14_2.seqs Order 0 >> output.test
./encoding-cost mark_f14_1.seqs mark_f14_2.seqs < counts_0_f_14_2 >> output.test
echo >> output.test
echo Train mark_f14_1.seqs Order 1 >>output.test
./encoding-cost mark_f14_2.seqs  mark_f14_1.seqs < counts_o1_f_14_1 >> output.test
echo >> output.test
echo Train mark_f14_2.seqs Order 1 >> output.test
./encoding-cost mark_f14_1.seqs  mark_f14_2.seqs < counts_o1_f_14_2 >> output.test
echo >> output.test
echo Constrained Alphabet Tests >> output.test
echo Train mark_f14_1.seqs Order 0 >> output.test
./encoding-cost --alphabet=ACDEFGHIKLMNPQRSTVWY mark_f14_2.seqs  mark_f14_1.seqs < counts_0_f_14_1_constrained >> output.test
echo >> output.test
echo Train mark_f14_2.seqs Order 0 >> output.test
./encoding-cost --alphabet=ACDEFGHIKLMNPQRSTVWY mark_f14_1.seqs  mark_f14_2.seqs < counts_0_f_14_2_constrained >> output.test
echo >> output.test
echo Train mark_f14_1.seqs Order 1 >>output.test
./encoding-cost --alphabet=ACDEFGHIKLMNPQRSTVWY mark_f14_2.seqs  mark_f14_1.seqs < counts_o1_f_14_1_constrained >> output.test
echo >> output.test
echo Train mark_f14_2.seqs Order 1 >> output.test
./encoding-cost --alphabet=ACDEFGHIKLMNPQRSTVWY mark_f14_1.seqs  mark_f14_2.seqs < counts_o1_f_14_2_constrained >> output.test
