# Stage 1: SAM (from samtools view -d lt) -> per-contig count of eligible ends with MAPQ>=MQ.
# Output: name <TAB> count, for every contig named by any lt record (count may be 0).
function bit(x,b){return int(x/b)%2}
function die(m){print "SUPPORTED_WINDOWS_FAILED: " m > "/dev/stderr"; bad=1; exit 1}
/^@/ || NF==0 {next}
{
  if (NF<11) die("Malformed SAM")
  n=0; for(i=12;i<=NF;i++) if (substr($i,1,5)=="lt:Z:") {n++; t=substr($i,6)}
  if (n==0) next
  if (n!=1) die("Duplicate lt tag")
  c=index(t,","); name=(c? substr(t,1,c-1) : t)
  if (!(name in q)) q[name]=0
  f=$2+0
  if (bit(f,4)||bit(f,256)||bit(f,512)||bit(f,1024)||bit(f,2048)) next
  if (bit(f,1) && !bit(f,2) && $7!="=" && $7!="*" && $7!=$3) next
  if (bit(f,1) && bit(f,64)==bit(f,128)) die("Bad segment flags")
  k=$1 SUBSEP (bit(f,64)*64+bit(f,128)*128)
  if (k in seen) die("Duplicate eligible end")
  seen[k]=1
  if ($5+0>=MQ) q[name]++
}
END{ if(bad) exit 1; for(x in q) print x "\t" q[x] }
