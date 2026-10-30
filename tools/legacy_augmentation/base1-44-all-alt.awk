# Select single-base REF/ALT records where every sample in columns 11..54 (samples 2..45)
# is homozygous ALT (1|1). Header and selected records are written to test.vcf.
BEGIN { OFS="\t" }
/^#/ { print > "test.vcf"; next }
{
    if(length($4)==1 && length($5)==1){
        ok=1
        for(i=11;i<=54;i++){
            split($i,a,":")
            if(a[1]!="1|1"){ ok=0; break }
        }
        if(ok) print > "test.vcf"
    }
}
