#!/bin/bash
set -euo pipefail
#conda activate /g/scb/zeller/karcher/anaconda3/envs/blast

# makeblastdb -in //g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_reviewed.fasta \
#   -dbtype prot \
#   -parse_seqids \
#   -title uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_reviewed \
#   -out /g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/blast_dbs/uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_reviewed

# makeblastdb -in //g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_all.fasta \
#   -dbtype prot \
#   -parse_seqids \
#   -title uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_all \
#   -out /g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/blast_dbs/uniprotkb_beta_glucuronidase_AND_taxono_2026_03_13_all

for strain_name in $(ls /g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/prokka)
    do
    for blast_db in $(ls /g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/blast_dbs | rev | cut -d "." -f2 | rev | sort | uniq )
    do
        echo "Processing $strain_name"
        echo "$blast_db"
        blastp -query /g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/prokka/${strain_name}/${strain_name}.faa \
        -db /g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/uniprot_data/blast_dbs/$blast_db \
        -outfmt 6 \
        -evalue 1e-5 \
        -num_threads 32 \
        -out /g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/blast_results/${strain_name}_${blast_db}_blast_results.tsv
        cat /g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/blast_results/${strain_name}_${blast_db}_blast_results.tsv | sed "s/^/${strain_name}___/" > a; 
        mv a /g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/revisions_2/blast_results/${strain_name}_${blast_db}_blast_results.tsv
    done
done