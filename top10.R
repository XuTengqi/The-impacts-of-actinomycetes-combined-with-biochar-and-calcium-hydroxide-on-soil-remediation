library(ggplot2)
a<-read.csv("E:/paperdata/top10.csv",row.names=1)
windowsFonts(HEL=windowsFont("Helvetica CE 55 Roman"),RMN=windowsFont("Times New Roman"),ARL=windowsFont("Arial"))
a$level<-factor(a$level,levels = c("CK","Act12","D74","B","Act12B","D74B","L","Act12L","D74L"))
a$AC<-factor(a$AC,levels = c("CK","B","L"))
a$A<-factor(a$A,levels = c("NI","Act12","D74"))
a$species<-factor(a$species,levels = c("Others","Nitrospirae","Planctomycetes","Gemmatimonadetes","Thermotogae","Chloroflexi","Bacteroidetes","Acidobacteria","Thaumarchaeota","Actinobacteria","Proteobacteria"))
colors<-c("black","sandybrown","cyan","seashell","yellow","grey", "blue","purple", "gold", "limegreen", "red")
p<-ggplot(a,aes(x=a$AC,y=a$value,fill=a$species))+geom_bar(position = "fill",stat = "identity")+theme_bw()+facet_grid(.~a$A)+scale_fill_manual(values=colors)+labs(x=" ", y="Relative Abundance",fill="Taxon")+
  + theme(text=element_text(family="RMN",size=14),axis.text.y=element_text(size=14),axis.text.x=element_text(size=14),legend.title=element_text(size=16),legend.text=element_text(size=16))+
  guides(fill=guide_legend(reverse=T,keywidth = 1, keyheight = 1))
p