#!/bin/sh

# smks="1_SnakeQC.smk 2_SnakeDemux.smk 3_SnakeMapping.smk 4_SnakeAssembly.smk 5_SnakeMismatch.smk 6_SnakeNascent.smk"
# smks="3_SnakeMapping.smk 4_SnakeAssembly.smk"

# mkdir -p logs
# for smk in $smks; do
#     echo "================================================================================"
#     echo $smk

#     if ( true ); then
#         log="logs/`basename $smk .smk`.log"
#         ./$smk > $log 2>&1
#     else
#         snakemake -s $smk -np
#     fi
# done
./3_SnakeMapping.smk
./6_SnakeMismatch.smk
./7_SnakeExpression.smk
