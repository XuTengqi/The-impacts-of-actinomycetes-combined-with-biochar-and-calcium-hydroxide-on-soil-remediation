#基于转化的RDA分析

library(vegan)
#读取数据
genus <- read.csv("F:/132/phylum_table.csv", row.names = 1)#这里我用了属水平绝对丰度
genust <- t(genus)#数据转置，使其和env格式一样
env <-  read.csv("F:/132/env_table.csv", row.names = 1)

#物种数据hellinger转化
genus_hel <- decostand(genust, method = 'hellinger')


#初步RDA分析

rda_tb.scaling2 <- summary(rda_tb, scaling = 2)


#R2校正
r2 <- RsquareAdj(rda_tb)
rda_adj <- r2$adj.r.squared 

#前向选择
vif.cca(rda_tb)#计算方差膨胀因子后发现有几个解释变量的VIF值较大，因此需要剔除部分变量
rda_tb_forward_r <- ordiR2step(rda(genus_hel~1, env, scale = FALSE), scope = formula(rda_tb), R2scope = rda_adj, direction = 'forward', permutations = 999)
plot(rda_tb_forward_r, scaling = 2, main = '前向选择后，II 型标尺', display = c('wa', 'cn'))#绘制双序图查看前向选择后的结果，“wa”代表使用物种加权和计算的样方坐标，“cn”代表约束成分（即解释变量）
summary(rda_tb_forward_r, scaling = 2)#查看前向选择后的summary统计信息（这里RDA1特征值占比0.60246，RDA2特征值占比0.24741）
RsquareAdj(rda_tb)$adj.r.squared#前向选择前R2校正结果
RsquareAdj(rda_tb_forward_r)$adj.r.squared#前向选择后R2校正结果（这里总解释量为0.4533907）
#RDA1解释量为0.60246*0.4533907*100%=27.31%，RDA2解释量为0.24741*0.4533907*100%=11.22%

#对约束轴进行置换检验
rda_tb_forward_r_test <- anova(rda_tb_forward_r, permutations = 999)
rda_tb_forward_r_test #整体
rda_tb_forward_r_test_axis <- anova(rda_tb_forward_r, by = 'axis', permutations = 999)
rda_tb_forward_r_test_axis  #逐个坐标轴
rda_tb_forward_r_test_axis$`Pr(>F)` <- p.adjust(rda_tb_forward_r_test_axis$`Pr(>F)`, method = 'bonferroni')
rda_tb_forward_r_test_axis #p值校正

#利用ggplot2作图

#提取样方和环境因子排序坐标，前两轴，II 型标尺
rda_tb_forward_r.scaling2 <- summary(rda_tb_forward_r, scaling = 2)
rda_tb_forward_r.site <- data.frame(rda_tb_forward_r.scaling2$sites)[1:2]
rda_tb_forward_r.env <- data.frame(rda_tb_forward_r.scaling2$biplot)[1:2]
rda_tb_forward_r.site$sample <- rownames(rda_tb_forward_r.site)
rda_tb_forward_r.env$sample <- rownames(rda_tb_forward_r.env)

#读取样本分组数据
group <- read.csv("D:/factor.csv", row.names = 1)

#合并样本分组信息，构建 ggplot2 作图数据集
rda_tb_forward_r.site <- merge(rda_tb_forward_r.site, group, by = "sample")


#作图
library(ggplot2)
rda_tb_forward_r.site$treatment<-factor(rda_tb_forward_r.site$treatment,levels=c("CK","B","F","L","A","AB","AF","AL"))#设置处理的顺序
rda_tb_forward_r.site$AC<-factor(rda_tb_forward_r.site$AC,levels=c("NI","A"))#设置处理的顺序
rda_tb_forward_r.site$AM<-factor(rda_tb_forward_r.site$AM,levels=c("CK","B","F","L"))#设置处理的顺序
col<-c('orange','red','green3','blue')#设置颜色
pch<-c(16,17)#设置图形形状

p <- ggplot(rda_tb_forward_r.site, aes(RDA1, RDA2)) + 
  geom_point( size = 3, aes(color = AM, shape = AC)) +
  scale_color_manual(values = col)+
  scale_shape_manual(values = pch)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = 'transparent')) +
  labs(size = 10, x = 'RDA1 (27.31%)', y = 'RDA2 (11.22%)') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = rda_tb_forward_r.env, aes(x = 0, y = 0, xend = RDA1,yend = RDA2), arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'black') +
  geom_text(data = rda_tb_forward_r.env, aes(RDA1 * 1.1, RDA2 * 1.1, label = sample), color = 'black', size = 3)

p
