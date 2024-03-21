## Tests for minimap2_index_modiifer
There are three tests to check base functionality.

### Common test
This test create two index files: regular and modified.
> minimap2 -d test/test.mni test/test.fasta  
> minimap2 -d test/test.modified.mni --vcf-file-with-variants test/test_long_chr1_not_the_same.vcf.gz test/test.fasta  
> diff test/test.mni test/test.modified.mni

Modified index file should contains extra string and some strings should be mismatched.

### Empty VCF test
If input VCF file contains no variants from reference modified index would be same as regular one.
> minimap2 -d test/test.mni test/test.fasta  
> minimap2 -d test/test.modified.mni --vcf-file-with-variants test/empty.vcf.gz test/test.fasta  
> diff test/test.mni test/test.modified.mni

Index files should be the same (empty output after third command).

### Pseudo-variants test
In this test `test/test_long_chr1_the_same.vcf.gz` contains pseudo-variants (like A -> A, C -> C, etc). This variants would be processed with no effect.
> minimap2 -d test/test.mni test/test.fasta  
> minimap2 -d test/test.modified.mni --vcf-file-with-variants test/test_long_chr1_the_same.vcf.gz test/test.fasta  
> diff test/test.mni test/test.modified.mni

Index files should be the same.
