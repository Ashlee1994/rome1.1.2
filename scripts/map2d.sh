source /apps/mpss/intel_parallel_2015/composer_xe_2015.1.133/bin/compilervars.sh intel64
source /apps/mpss/intel_parallel_2015/impi/5.0.2.044/bin64/mpivars.sh

input_fn="./Inf_data1"
output_fn="./Class2D/Inf_data1_ml2d/Inf_data1_m2d_K100"
map2d_classes=100
offset_range=10
offset_step=2
psi_step=10
map2d_iter=30
pixel_size=1.72
nr_pool=20
others="-ctf -norm"
mpirun -n 32 -f ./all_machines_ib -perhost 1  ./bin/rome_map -i $input_fn -o $output_fn -K $map2d_classes -angpix $pixel_size -iter $map2d_iter -pool $nr_pool -offset_range $offset_range -offset_step $offset_step -psi_step $psi_step $others > map2d_K100_Inf_data1.r
