source /apps/mpss/intel_parallel_2015/composer_xe_2015.1.133/bin/compilervars.sh intel64
source /apps/mpss/intel_parallel_2015/impi/5.0.2.044/bin64/mpivars.sh
input_fn="./RP_data1"
output_fn="./Class2D/RP_data1_sml/RP_data1_sml_K50"
nr_classes=50
pixel_size=2.00
nr_iter=30
others="-ctf"
mpirun -n 32 -f ./all_machines_ib -perhost 1  ./bin/rome_sml -i $input_fn -o $output_fn -K $nr_classes -angpix $pixel_size -iter $nr_iter $others
