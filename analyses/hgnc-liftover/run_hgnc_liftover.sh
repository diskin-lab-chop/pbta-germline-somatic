# Run HGNC liftover on CPG list

# define scripts dir 
scripts_dir=(../../scripts/)

# define input files 
cpgs=(input/cpg.txt)
complete=(input/hgnc_complete_set_2023-04-01.txt)

# Define output file
out_cpg=(results/cpg_liftover.tsv)

# Run liftover
python $scripts_dir/update_fusion_gene_symbols.py -g $complete -f $cpgs -o $out_cpg -u "Gene"

# extract CPG name column
awk -F "\t" '{print $1}' $out_cpg  | sed '1d' > results/cpg.txt
