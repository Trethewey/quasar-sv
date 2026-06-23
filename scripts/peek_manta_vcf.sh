#!/usr/bin/env bash
# Quick inspector: print Manta VCF records hitting BCL6, IGH, IGL, IGK, MYC, BCL2 windows.
VCF=$1
# Auto-detect gzip vs plain text
case "$VCF" in *.gz) CAT="zcat" ;; *) CAT="cat" ;; esac
$CAT "$VCF" | awk 'BEGIN{OFS="\t"}
/^#/ { next }
{
  chr=$1; pos=$2; vid=$3; alt=$5;
  if ((chr=="chr3"  && pos>=187000000 && pos<=188000000) ||
      (chr=="chr14" && pos>=105500000 && pos<=107000000) ||
      (chr=="chr22" && pos>=22000000  && pos<=23000000)  ||
      (chr=="chr2"  && pos>=88800000  && pos<=90300000)  ||
      (chr=="chr8"  && pos>=127500000 && pos<=128100000) ||
      (chr=="chr18" && pos>=63100000  && pos<=63400000)) {
    print chr, pos, vid, substr(alt, 1, 60), $7
  }
}'
