
cp /var/tmp/XSBench_aocx/base/CalculateXS.aocx .
./XSBench -t 15000000 | tee results/base.txt
start_line=$(grep 'OCL_KERNEL_TS_SEC' results/base.txt)
end_line=$(grep 'OCL_POST_TS_SEC' results/base.txt)
start_time=`echo "$start_line" | sed 's/^.*OCL_KERNEL_TS_SEC=//'`
end_time=`echo "$end_line" | sed 's/^.*OCL_POST_TS_SEC=//'`
ssh knight ~/gitwork/ocl-iabench/wt310script/wt310samples_sqlite.py $start_time $end_time ~/base.txt 

cp /var/tmp/XSBench_aocx/cc2mgridpoint/CalculateXS.aocx .
./XSBench -t 15000000 | tee results/cc2mgridpoint.txt
start_line=$(grep 'OCL_KERNEL_TS_SEC' results/cc2mgridpoint.txt)
end_line=$(grep 'OCL_POST_TS_SEC' results/cc2mgridpoint.txt)
start_time=`echo "$start_line" | sed 's/^.*OCL_KERNEL_TS_SEC=//'`
end_time=`echo "$end_line" | sed 's/^.*OCL_POST_TS_SEC=//'`
ssh knight ~/gitwork/ocl-iabench/wt310script/wt310samples_sqlite.py $start_time $end_time ~/cc2mgridpoint.txt

type=fma 
cp /var/tmp/XSBench_aocx/$type/CalculateXS.aocx .
./XSBench -t 15000000 | tee results/$type.txt
start_line=$(grep 'OCL_KERNEL_TS_SEC' results/$type.txt)
end_line=$(grep 'OCL_POST_TS_SEC' results/$type.txt)
start_time=`echo "$start_line" | sed 's/^.*OCL_KERNEL_TS_SEC=//'`
end_time=`echo "$end_line" | sed 's/^.*OCL_POST_TS_SEC=//'`
ssh knight ~/gitwork/ocl-iabench/wt310script/wt310samples_sqlite.py $start_time $end_time ~/$type.txt

type=optbsearch
cp /var/tmp/XSBench_aocx/$type/optbCalculateXS.aocx .
./XSBench_optb -t 15000000 | tee results/$type.txt
start_line=$(grep 'OCL_KERNEL_TS_SEC' results/$type.txt)
end_line=$(grep 'OCL_POST_TS_SEC' results/$type.txt)
start_time=`echo "$start_line" | sed 's/^.*OCL_KERNEL_TS_SEC=//'`
end_time=`echo "$end_line" | sed 's/^.*OCL_POST_TS_SEC=//'`
ssh knight ~/gitwork/ocl-iabench/wt310script/wt310samples_sqlite.py $start_time $end_time ~/$type.txt

type=vectornucs
cp /var/tmp/XSBench_aocx/$type/optbCalculateXS.aocx .
./XSBench_optb_vector -t 15000000 | tee results/$type.txt
start_line=$(grep 'OCL_KERNEL_TS_SEC' results/$type.txt)
end_line=$(grep 'OCL_POST_TS_SEC' results/$type.txt)
start_time=`echo "$start_line" | sed 's/^.*OCL_KERNEL_TS_SEC=//'`
end_time=`echo "$end_line" | sed 's/^.*OCL_POST_TS_SEC=//'`
ssh knight ~/gitwork/ocl-iabench/wt310script/wt310samples_sqlite.py $start_time $end_time ~/$type.txt
