library(ggplot2)
windowsFonts(HEL=windowsFont("Helvetica CE 55 Roman"),RMN=windowsFont("Times New Roman"),ARL=windowsFont("Arial"))
a<-read.csv("E:/paperdata/Cdform.csv",row.names = 1)
a$level<-factor(a$level,levels = c("CK","Act12","D74","B","Act12B","D74B","L","Act12L","D74L"))
a$Ac<-factor(a$Ac,levels = c("NI","Act12","D74"))
a$Am<-factor(a$Am,levels = c("CK","B","L"))
a$shape<-factor(a$shape,levels = c("Residual form","Organic bound form","Fe-Mn oxide bound form","Carbonate bound form","Exchangeable form"))
colors<-c("blue","purple", "gold", "limegreen", "red")
p<-ggplot(a, aes(x = Ac, y = value, fill = shape))+geom_bar(position = "fill", stat = "identity")
g<-p+theme_bw()+facet_grid(.~Am)+scale_fill_manual(values=colors)+labs(x=" ", y="Relative Abundance",fill="Taxon")+
g + theme(text=element_text(family="RMN",size=14),axis.text.y=element_text(size=14),axis.text.x=element_text(size=14),legend.title=element_text(size=16),legend.text=element_text(size=16))+
  guides(fill=guide_legend(reverse=T,keywidth = 1, keyheight = 1))