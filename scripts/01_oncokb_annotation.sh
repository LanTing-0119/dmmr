# https://oncotree.mskcc.org/#/home
IMAF38="/Users/lanting/Project/dmmr/data/vep_merge.maf"
OMAF38="/Users/lanting/Project/dmmr/data/vep_merge_oncokb.maf"
IC="/Users/lanting/Project/dmmr/outputs/oncokb_clinical.csv"
TOKEN="84a75c6f-7daf-4fee-86ce-ae3f12f31413" #OncoKB API Token
MafAnnotator=/Users/lanting/Project/dmmr/scripts/tools/oncokb-annotator-master/MafAnnotator.py
conda activate oncokb_annotation
cd /Users/lanting/Project/dmmr/scripts/tools/oncokb-annotator-master/
python3 ${MafAnnotator} -i "$IMAF38" -o "$OMAF38" -c "$IC" -b "$TOKEN" -t 'COADREAD'
# without pointing out cancer type will affect Therapeutics/Diagnostic/Prognostic implications but not
# python3 ${MafAnnotator} -i "$IMAF38" -o "$OMAF38" -c "$IC" -b "$TOKEN"
