# Stage 3: events sorted by (chrom_index, pos, dl) -> merged windows [pos-W, pos+dl+W) clipped to the contig,
# clusters longer than MAXLEN split recursively at the widest gap between event intervals (else midpoint).
function piece(s,e,   i,a,b,best,cut,reach,have){
  if (e-s<=MAXLEN) { printf "%s\t%d\t%d\n", CH, s, e; return }
  best=0; have=0; cut=-1
  for(i=1;i<=n;i++){
    if (EA[i]<e && EB[i]>s) {
      a=(EA[i]>s?EA[i]:s); b=(EB[i]<e?EB[i]:e)
      if (have && a-reach>best) { best=a-reach; cut=int((reach+a)/2) }
      if (!have || b>reach) reach=b
      have=1
    }
  }
  if (cut<0 || cut<=s || cut>=e) cut=int((s+e)/2)
  piece(s,cut); piece(cut,e)
}
function flush(){ if (n) { piece(CS,CE) }; n=0 }
BEGIN{ while((getline l < FAI)>0){ split(l,f,"\t"); L[f[1]]=f[2]+0 } }
{
  ch=$4; pos=$2+0; dl=$3+0
  s=pos-W; if (s<0) s=0; e=pos+dl+W; if (e>L[ch]) e=L[ch]
  if (n && ch==CH && s<=CE) { if (e>CE) CE=e }
  else { flush(); CH=ch; CS=s; CE=e }
  n++; EA[n]=pos-1; EB[n]=pos+dl
}
END{ flush() }
