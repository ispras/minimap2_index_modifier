# Select single-base REF/ALT records with INFO AF > 0.75.
# Header and selected records are written to test.vcf.
BEGIN { OFS="\t" }
/^#/ { print > "test.vcf"; next }
{
    if (length($4) == 1 && length($5) == 1) {
        af_val = ""
        n = split($8, info, ";")
        for (i=1; i<=n; i++) {
            if (info[i] ~ /^AF=/) {
                af_val = substr(info[i], 4)
                break
            }
        }
        if (af_val != "" && af_val > 0.75) {
            print > "test.vcf"
        }
    }
}
