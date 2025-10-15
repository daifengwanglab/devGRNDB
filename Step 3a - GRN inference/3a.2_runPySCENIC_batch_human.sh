 
 
 #!/bin/bash

species=("Human")
cases=("Forebrain")
lineages=("Broad")

# Loop by index
for i in "${!cases[@]}"; do
  sp=${species[$i]}
  case=${cases[$i]}
  lin=${lineages[$i]}

  echo "Running pySCENIC for $sp - $case - $lin"

  # Step 1: GRN
  pyscenic grn \
    --method grnboost2 \
    "${sp,,}/${case}/${case}_${lin}.loom" \
    allTFs_hg38.txt \
    -o "${sp,,}/${case}/${case}_${lin}_adj.csv" \
    --num_workers 30

  # Step 2: ctx
  pyscenic ctx \
    "${sp,,}/${case}/${case}_${lin}_adj.csv" \
    hg38v10pyscenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
    hg38v10pyscenic/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
    hg38v10pyscenic/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname hg38v10pyscenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname "${sp,,}/${case}/${case}_${lin}.loom" \
    --output "${sp,,}/${case}/${case}_${lin}_ctx.csv" \
    --num_workers 30

  echo "Done: $sp - $case - $lin"
done
