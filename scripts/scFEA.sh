#! /bin/bash

#input file
INFILE=$1
#output directory
OUTDIR=$2

# confirm right number of arguments
if [ $# -le 1 -o $# -ge 3 ]; then
	echo "Wrong number of arguments. Provide 2 arguments [input file] [output directory]"
	exit 1
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $SCRIPT_DIR

python scFEA.py --data_dir data --input_dir .. \
                    --test_file $INFILE \
                    --res_dir ../$OUTDIR \
                    --moduleGene_file module_gene_complete_mouse_m168.csv \
                    --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv \
                    --sc_imputation True \
                    --output_flux_file ../data/Flux.csv \
                    --output_balance_file ../data/Metabolome_balance.csv
