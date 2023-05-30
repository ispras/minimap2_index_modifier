## Build
install htslib 1.17
make


## Run modified index creation
bgzip -c test.vcf > test.vcf.gz
tabix -p vcf test.vcf.gz
./minimap2 -d test/test.modified.mni --vcf-file-with-variants test/test.vcf.gz test/test.fasta
