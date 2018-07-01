source /apps/mpss/intel_parallel_2015/composer_xe_2015.1.133/bin/compilervars.sh intel64
source /apps/mpss/intel_parallel_2015/impi/5.0.2.044/bin64/mpivars.sh

input_fn="./Inf_data2"
output_fn="./Class2D/Inf_data2_deep2d/Inf_data2_K300"
deep2d_classes=300
map2d_classes=50
offset_range=10
offset_step=2
psi_step=10
map2d_iter=30
sml_iter=30
pixel_size=1.72
nr_pool=50
others="-ctf -norm"
mpirun -n 32 -f ./all_machines_ib -perhost 1  ./bin/rome_deep2d -i $input_fn -o $output_fn -map2d_K $map2d_classes -sml_K $deep2d_classes -angpix $pixel_size -map2d_iter $map2d_iter -sml_iter $sml_iter -pool $nr_pool -offset_range $offset_range -offset_step $offset_step -psi_step $psi_step $others > Inf_data2_K300.r
