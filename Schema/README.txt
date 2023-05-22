This folder contains the Txo2 (6J1P) pdb file, parent aligned file (Pa11_6P1J_C1-A10.aln) and multiple sequence alignment file (Pa11_MSA_C1-A10.aln) of the second module of PvdD to the 8 alternate CA sequences. 

A contacts file was created by running the command:
python schema-tools\schemacontacts.py -pdb 6p1j.pdb -msa Pa11_MSA_C1-A10.aln -pdbal Pa11_6P1J_C1-A10.aln -o Pa11_contacts.txt

The contacts and MSA are then used by Schema_profile_Pvdd_motifs.py to produce a SCHEMA plot of the 8 sequences (Profile_Pa11_Txo2_labelled.png), labelled by substrate specificity. Conserved motifs, sites 1, 2 and X; are are labelled by dashed or solid lines.