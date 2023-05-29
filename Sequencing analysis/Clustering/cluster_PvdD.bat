usearch11.0.667_win32.exe -fastx_uniques processed_reverse_sequences.fasta -fastaout processed_reverse_sequences_uniques.fasta
usearch11.0.667_win32.exe -sortbylength processed_reverse_sequences_uniques.fasta -fastaout processed_reverse_sequences_sorted.fasta
usearch11.0.667_win32.exe -cluster_fast processed_reverse_sequences_sorted.fasta -id 0.97 -centroids processed_reverse_sequences_otus.fasta -uc processed_reverse_sequences_renamed.uc -relabel OTU
usearch11.0.667_win32.exe -otutab processed_reverse_sequences.fasta -otus processed_reverse_sequences_otus.fasta -id 0.95 -otutabout processed_reverse_sequences_otutab.txt -mapout processed_reverse_sequences_map.txt
PAUSE