
#python2 getCodonAversionChanges.py -i ${1}  -b -a clinvar_codon_aversion -o clinvar_codon_pairing
python2 scripts/getCodonAversionChanges.py -i ${1}  -b -a ${1}_codon_aversion -o ${1}_codon_pairing
