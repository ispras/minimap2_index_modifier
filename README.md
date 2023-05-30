## Build
install htslib 1.17
make

## Make gz and gz.tbi from VCF
bgzip -c filename.vcf > filename.vcf.gz
tabix -p vcf filename.vcf.gz

## Run modified index creation
bgzip -c test.vcf > test.vcf.gz
tabix -p vcf test.vcf.gz
./minimap2 -d test/test.modified.mni --vcf-file-with-variants test/test_long_chr1.vcf.gz test/test.fasta

## Run test

### Empty test
./minimap2 -d test/test.mni test/test.fasta
./minimap2 -d test/test.modified.mni --vcf-file-with-variants test/empty.vcf.gz test/test.fasta
diff test/test.mni test/test.modified.mni


### The same test
./minimap2 -d test/test.mni test/test.fasta
./minimap2 -d test/test.modified.mni --vcf-file-with-variants test/test_long_chr1_the_same.vcf.gz test/test.fasta
diff test/test.mni test/test.modified.mni

### Not the same test
./minimap2 -d test/test.mni test/test.fasta
./minimap2 -d test/test.modified.mni --vcf-file-with-variants test/test_long_chr1_not_the_same.vcf.gz test/test.fasta
diff test/test.mni test/test.modified.mni
