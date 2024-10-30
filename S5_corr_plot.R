library(ggplot2)
library(pheatmap)
library(corrplot)

mode='count'
fr=read.table('DGE_count.txt',header=T,row.names=1,sep='\t')
wl=read.table('whole.DGE_count.txt',header=T,row.names=1,sep='\t')
e=cbind(fr,wl)

aver.data=data.frame(Fraction=rowMeans(fr),Whole=rowMeans(wl))
pdata=apply(aver.data,2,function(x){log(x+0.25)/log(2)})
pdata=data.frame(pdata)
aver.corr=cor(pdata$Fraction,pdata$Whole)
aver.celln=dim(pdata)[1]
corr_celln=data.frame(corr=aver.corr,cellnum=aver.celln)
write.table(t(corr_celln),'F-W_corr.txt',col.names=F)
p<-ggplot(pdata, aes(x=Fraction,y=Whole))+geom_point()+geom_smooth(color="red",method="lm",se=F)
p<-p+theme_bw()+
   geom_text(label=paste('R=',round(aver.corr,3)),x=max(pdata$Fraction)*0.08,y=max(pdata$Whole)*0.9,size=5)+
   labs(x="Log2(count+0.25) of scBiopsy-seq", y="Log2(count+0.25) of scRNA-seq")+
   theme(panel.border=element_blank(),axis.line=element_line(color="black"),
   panel.grid.major=element_blank(),panel.grid.minor=element_blank())
pdf('F-W_corr.pdf',width=4,height=4)
print(p)
dev.off()
