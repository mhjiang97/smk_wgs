#!/usr/bin/env bash
# shellcheck disable=SC2154

set -x

{ vcf=${snakemake_input[vcf]}
vcf_filtered=${snakemake_output[vcf]}
min_size=${snakemake_params[min_size]}
min_reads=${snakemake_params[min_reads]}
min_coverage=${snakemake_params[min_coverage]}
min_dhffc=${snakemake_params[min_dhffc]}
max_dhbfc=${snakemake_params[max_dhbfc]}
caller=${snakemake_wildcards[caller_sv]}

formula_length="INFO/SVTYPE != \"BND\" & ABS(INFO/SVLEN) < ${min_size}"
formula_del="INFO/SVTYPE == \"DEL\" & FMT/DHFFC[0] > ${min_dhffc}"
formula_dup="INFO/SVTYPE == \"DUP\" & FMT/DHBFC[0] < ${max_dhbfc}"

if [ "${caller}" == "gridss" ]; then

    formula_reads="((INFO/ASRP + INFO/ASSR) < ${min_reads}) & ((INFO/SR + INFO/RP) < ${min_reads}) & ((INFO/SR + INFO/BSC) < ${min_reads})"
    formula_coverage="((INFO/ASSR + INFO/REF + INFO/ASRP + INFO/REFPAIR) < ${min_coverage}) & ((INFO/SR + INFO/RP + INFO/REF + INFO/REFPAIR) < ${min_coverage}) & ((INFO/SR + INFO/BSC + INFO/REF) < ${min_coverage})"

elif [ "${caller}" == "manta" ]; then

    formula_reads="FMT/PR[0:1] + FMT/SR[0:1] < ${min_reads}"
    formula_coverage="((SUM(FMT/PR[0:]) + SUM(FMT/SR[0:])) < ${min_coverage}) | (INFO/SVTYPE == \"BND\" & (INFO/BND_DEPTH < ${min_coverage}) | (INFO/MATE_BND_DEPTH < ${min_coverage}))"

elif [ "${caller}" == "svaba" ]; then

    formula_reads="FMT/AD[0] < ${min_reads}"
    formula_coverage="FMT/DP[0] < ${min_coverage}"

elif [ "${caller}" == "tiddit" ]; then

    formula_reads="(FMT/DV[0] + FMT/RV[0]) < ${min_reads}"
    formula_coverage="(FMT/DV[0] + FMT/RV[0] + SUM(FMT/RR[0:])) < ${min_coverage}"

elif [ "${caller}" == "wham" ]; then

    formula_reads="INFO/A < ${min_reads} | FMT/SP[0] < ${min_reads}"
    formula_coverage="FMT/SP[0] < 0"

else
    echo "$(date +"%Y-%m-%d %H:%M:%S") [ERROR] Unknown caller: ${caller}"
    exit 1
fi

formula="(${formula_length}) | (${formula_coverage}) | (${formula_reads}) | (${formula_del}) | (${formula_dup})"

grep -E "^#|^chr" "${vcf}" \
    | bcftools filter -e "${formula}" -Ov - \
    > "${vcf_filtered}"; } \
1> "${snakemake_log[0]}" 2>&1
