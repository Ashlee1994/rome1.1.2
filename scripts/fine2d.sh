source /apps/mpss/intel_parallel_2015/composer_xe_2015.1.133/bin/compilervars.sh intel64
source /apps/mpss/intel_parallel_2015/impi/5.0.2.044/bin64/mpivars.sh
input_fn="./deep2d_Inf_data2.star"
output_fn="./Class2D/Inf_data2_deep2d_fine/"
pixel_size=1.72
nr_iter=30
mpirun -n 32 -f ./all_machines_ib -perhost 1  ./ragon05/yb/rome10a_release/bin/rome_sml -search class.config -i $input_fn -o $output_fn -angpix $pixel_size -iter $nr_iter > result.r
#ring_all_top_side_fine2d.r
