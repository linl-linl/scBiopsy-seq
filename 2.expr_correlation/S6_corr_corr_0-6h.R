library(ggplot2)
library(pheatmap)
library(corrplot)

dge=read.table('DGE_count.txt',header=T,row.names=1,sep='\t')
anno=read.table('anno.txt',header=T,row.names=1)

t0=rownames(anno)[anno$TIME=='t0']
t1=rownames(anno)[anno$TIME=='t1']
t4=rownames(anno)[anno$TIME=='t4']
t6=rownames(anno)[anno$TIME=='t6']
aver.data=data.frame(t.0h=rowMeans(dge[,t0]), t.1h=rowMeans(dge[,t1]),
                     t.4h=rowMeans(dge[,t4]), t.6h=rowMeans(dge[,t6]))
e=aver.data

corr=e[1:dim(e)[2],]
rownames(corr)=colnames(corr)
ngene=corr
for (i in colnames(e)){
  for (j in colnames(e)){
    ei=log(e[,i]+0.25)/log(2)
    ej=log(e[,j]+0.25)/log(2)
    corr[i,j]=cor(ei,ej)
    ngene[i,j]=length(ei)
    if (which(colnames(e)==i)<which(colnames(e)==j)){
      pdata=data.frame(ei=ei,ej=ej)
      p<-ggplot(pdata, aes(x=ei,y=ej))+geom_point()+geom_smooth(color="red",method="lm",se=F)
      p<-p+theme_bw()+
           geom_text(label=paste('R=',round(corr[i,j],3)),x=max(ei)*0.08,y=max(ej)*0.9,size=5)+
           labs(x=paste0("Log2(count+0.25) of ",i), y=paste0("Log2(count+0.25) of ",j))+
           theme(panel.border=element_blank(),axis.line=element_line(color="black"),
           panel.grid.major=element_blank(),panel.grid.minor=element_blank())
      pdf(paste0('aver_corr.',gsub('t.','',i),'-',gsub('t.','',j),'.pdf'),width=4,height=4)
      print(p)
      dev.off()
    }
  }
}
write.table(corr,'aver_corr.txt',sep='\t')
pdf('aver_corr.pdf',width=6,height=6)
p.cor<-corrplot(as.matrix(corr),method='color',type='full',order='AOE',
         cl.pos='r',cl.ratio=0.2,cl.cex=0.8,cl.align.text='l',
         tl.pos='lt',tl.cex=1,tl.col='black',
         addCoef.col='white',number.cex=0.9,mar=c(1,1,1,0))
dev.off()
write.table(ngene,'aver_corr.nGene.txt',sep='\t')
pdf('aver_corr.nGene.pdf',width=3,height=3)
sample_order=unique(p.cor$corrPos$xName)
pheatmap(corr[sample_order,sample_order],display_numbers=ngene[sample_order,sample_order],
         color='white',number_color='black',fontsize_number=8,
         cluster_rows=F,cluster_cols=F,legend=F,
         cellwidth=30,cellheight=30,border_color='black',fontsize_row=12,fontsize_col=12,angle_col=90)
dev.off()
