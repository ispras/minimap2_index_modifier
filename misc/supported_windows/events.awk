# Stage 2: FAI + support counts + panel JSONL -> deduplicated events of contigs with >=N supporting ends.
# Output: chrom_index <TAB> pos <TAB> len(ref)-1 <TAB> chrom <TAB> ref <TAB> alt  (chrom_index = FAI order; pos = 1-based POS)
# Validates the whole panel (unless -v VALIDATE=0: only selected templates); expects canonical json.dumps lines.
function die(m){print "SUPPORTED_WINDOWS_FAILED: " m > "/dev/stderr"; bad=1; exit 1}
function field(key,  s){ s=index($0,"\"" key "\": "); if(!s) die("Panel line not canonical"); return substr($0,s+length(key)+4) }
FILENAME==FAI {L[$1]=$2+0; ord[$1]=++nf; next}
FILENAME==SUP {q[$1]=$2+0; next}
{
  if (substr($0,1,10)!="{\"name\": \"") die("Panel line not canonical")
  name=substr($0,11); name=substr(name,1,index(name,"\"")-1)
  if (name in dup) die("Duplicate template name"); dup[name]=1
  if (VALIDATE=="0" && !((name in q) && q[name]>=N)) next
  ch=field("chrom"); ch=substr(ch,2,index(substr(ch,2),"\"")-1)
  o=field("origins"); o=substr(o,2,index(o,"]")-2)
  if (o=="") die("Unknown carrier")
  no=split(o,oa,", "); for(i=1;i<=no;i++) if (oa[i]!="\"hap1\"" && oa[i]!="\"hap2\"") die("Unknown carrier")
  ev=field("events"); ev=substr(ev,2,index(ev,"]")-2)
  if (ev=="") die("Empty template")
  sel=(name in q) && q[name]>=N
  ne=split(ev,ea,"}, {")
  for(i=1;i<=ne;i++){
    e=ea[i]; gsub(/[{}]/,"",e)
    if (split(e,p,"\"")!=15 || p[2]!="chrom" || p[6]!="pos" || p[8]!="ref" || p[12]!="alt") die("Panel event not canonical")
    ech=p[4]; pos=p[7]; gsub(/[:, ]/,"",pos); ref=p[10]; alt=p[14]
    if (pos !~ /^-?[0-9]+$/) die("Panel event not canonical")
    pos+=0
    if (ech!=ch || !(ech in L) || ref=="" || alt=="" || pos<1 || pos>L[ech] || pos-1+length(ref)>L[ech] || length(ref)==length(alt)) die("Bad event geometry")
    if (sel) { k=ech SUBSEP pos SUBSEP ref SUBSEP alt; if (!(k in out)) { out[k]=1; print ord[ech] "\t" pos "\t" length(ref)-1 "\t" ech "\t" ref "\t" alt } }
  }
}
END{ if(bad) exit 1; for(x in q) if(!(x in dup)) die("Unknown lt") }
