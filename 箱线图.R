library(reshape2)
alpha <- read.csv("F:/132/alpha.csv",stringsAsFactors = FALSE)
alpha$group2 <- factor(alpha$group2)
alpha1 <- melt(alpha, id = c('samples', 'group1', 'group2'))
alpha2 <- subset(alpha1, variable == 'chao1')
alpha3 <- subset(alpha2, group1 == 'c')
library(ggplot2)

p <- ggplot(alpha1, aes(x = group2, y = value, fill = group1)) + 
  geom_boxplot(outlier.size = 0.3, size = 0.25) +
  facet_wrap(~variable, 2, scales = 'free') +
  labs(x = '', y = '') +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())
      axis.alpha.y=element_text(size=14)
      axis.alpha.x=element_text(size = 16)
      axis.ticks.x=element_blank( )
      legend.title=element_text(face="plain",size=14)
      legend.text=element_text(size=12)

#ggsave('ggplot2.box_facet.pdf', p, width = 7, height = 5)
ggsave('ggplot2.box_facet.png', p, width = 7, height = 5)

p

