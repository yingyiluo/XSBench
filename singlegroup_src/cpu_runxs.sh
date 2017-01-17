type=base
spec=Intel1-protos
cp kernel/$type.cl CalculateXS.cl
numactl -m 1 --cpunodebind=1 ./XSBench | tee cpu_results/$type-$spec.txt

type=cc
spec=Intel1-protos
cp kernel/$type.cl CalculateXS.cl
numactl -m 1 --cpunodebind=1 ./XSBench | tee cpu_results/$type-$spec.txt

type=fma
spec=Intel1-protos
cp kernel/$type.cl CalculateXS.cl
numactl -m 1 --cpunodebind=1 ./XSBench | tee cpu_results/$type-$spec.txt

type=optb
spec=Intel1-protos
cp kernel/$type.cl CalculateXS.cl
numactl -m 1 --cpunodebind=1 ./XSBench_optb | tee cpu_results/$type-$spec.txt

type=vectoroptb
spec=Intel1-protos
cp kernel/$type.cl CalculateXS.cl
numactl -m 1 --cpunodebind=1 ./XSBench_optb_vector | tee cpu_results/$type-$spec.txt

