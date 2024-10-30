
system("ls *_count.txt > list_cnt")
system("ls *_fpkm.txt > list_fpkm")
mt = read.table('/home/disk/linl/reference/GRCh38/mt_id.txt')
mt = as.character(mt$V1)
rp = read.table('/home/disk/linl/dmf_live_seq/live_seq/ncbi_rp_id.txt')
rp = as.character(rp$V1)

merge_dge <- function(list_expr, dge_name){
  expr = read.table(list_expr)
  expr = as.character(expr$V1)
  n = 0
  for (i in expr){
    n = n+1
    e = read.table(i, header=T, row.names=1, sep='\t')
    if (n==1){
      dge = e
    }else{
      dge = cbind(dge,e)
    }
  }
  write.table(dge, dge_name, sep='\t')
  dge_clean = dge[which( !(rownames(dge) %in% mt) ), ]
  dge_clean = dge_clean[which( !(rownames(dge_clean) %in% rp) ), ]
  dge_clean_name = gsub('.txt', '_clean.txt', dge_name)
  write.table(dge_clean, dge_clean_name, sep='\t')
}

merge_dge('list_cnt', 'DGE_count.txt')
merge_dge('list_fpkm', 'DGE_fpkm.txt')
system("rm list_cnt")
system("rm list_fpkm")

