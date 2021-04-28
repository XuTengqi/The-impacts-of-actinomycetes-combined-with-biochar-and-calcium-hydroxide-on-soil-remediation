library(vegan)

##读取数据
#读入物种数据，以细菌 OTU 水平丰度表为例
otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))

#读取环境数据
env <- read.delim('env_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#############################
##根据原理一步步计算 db-RDA
#计算样方距离，以 Bray-curtis 距离为例，详情 ?vegdist
dis_bray <- vegdist(otu, method = 'bray')

#或者直接使用现有的距离矩阵，这里同样为 Bray-curtis 距离
dis_bray <- as.dist(read.delim('bray_distance.txt', row.names = 1, sep = '\t', check.names = FALSE))

#PCoA 排序，这里通过 add = TRUE校正负特征值，详情 ?cmdscale
pcoa <- cmdscale(dis_bray, k = nrow(otu) - 1, eig = TRUE, add = TRUE)

#提取 PCoA 样方得分（坐标）
pcoa_site <- pcoa$point

#db-RDA，环境变量与 PCoA 轴的多元回归
#通过 vegan 包的 RDA 函数 rda() 执行，详情 ?rda
db_rda <- rda(pcoa_site, env, scale = FALSE)

#被动拟合物种得分
v.eig <- t(otu) %*% db_rda$CCA$u/sqrt(nrow(otu) - 1)
db_rda$CCA$v <- decostand(v.eig, 'normalize', MARGIN = 2)
v.eig <- t(otu) %*% db_rda$CA$u/sqrt(nrow(otu) - 1)
db_rda$CA$v <- decostand(v.eig, 'normalize', MARGIN = 2)

#先作图展示下，详情 ?plot.cca
#样方展示为点，物种暂且展示为“+”，环境变量为向量
par(mfrow = c(1, 2))

plot(db_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1, main = 'I型标尺，双序图')
points(db_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), cex = 1)
text(db_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'blue', cex = 0.8)

plot(db_rda, type = 'n', display = c('wa', 'sp', 'cn'), choices = 1:2, scaling = 1, main = 'I型标尺，三序图')
points(db_rda, choices = 1:2, scaling = 1, display = 'sp', pch = 3, col = 'gray', cex = 1)
points(db_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), cex = 1)
text(db_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'blue', cex = 0.8)

#############################
#或者，capscale() 提供了直接运行的方法，详情 ?capscale

#输入原始数据，指定距离类型，例如样方距离以 Bray-curtis 距离为例
#计算结果中包含物种得分
db_rda <- capscale(otu~., env, distance = 'bray', add = TRUE)

#或者直接输入已经计算/或读入的距离测度
db_rda <- capscale(dis_bray~., env, add = TRUE)
#但是这种情况下，物种得分丢失，需要被动添加
v.eig <- t(otu) %*% db_rda$CCA$u/sqrt(nrow(otu) - 1)
db_rda$CCA$v <- decostand(v.eig, 'normalize', MARGIN = 2)
v.eig <- t(otu) %*% db_rda$CA$u/sqrt(nrow(otu) - 1)
db_rda$CA$v <- decostand(v.eig, 'normalize', MARGIN = 2)

#先作图展示下
#样方展示为点，物种暂且展示为“+”，环境变量为向量
par(mfrow = c(1, 2))

plot(db_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1, main = 'I型标尺，双序图')
points(db_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), cex = 1)
text(db_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'blue', cex = 0.8)

plot(db_rda, type = 'n', display = c('wa', 'sp', 'cn'), choices = 1:2, scaling = 1, main = 'I型标尺，三序图')
points(db_rda, choices = 1:2, scaling = 1, display = 'sp', pch = 3, col = 'gray', cex = 1)
points(db_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), cex = 1)
text(db_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'blue', cex = 0.8)

#############################
#以及 dbrda()，详情 ?dbrda

#输入原始数据，指定距离类型，例如样方距离以 Bray-curtis 距离为例
db_rda <- dbrda(otu~., env, distance = 'bray', add = TRUE)

#或者直接输入已经计算/或读入的距离测度
db_rda <- dbrda(dis_bray~., env, add = TRUE)

#两种情况下，dbrda() 均默认不计算物种得分，可通过 sppscores() 添加，详情 ?sppscores
sppscores(db_rda) <- otu

#先作图展示下
#样方展示为点，物种暂且展示为“+”，环境变量为向量
par(mfrow = c(1, 2))

plot(db_rda, type = 'n', display = c('wa', 'cn'), choices = 1:2, scaling = 1, main = 'I型标尺，双序图')
points(db_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), cex = 1)
text(db_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'blue', cex = 0.8)

plot(db_rda, type = 'n', display = c('wa', 'sp', 'cn'), choices = 1:2, scaling = 1, main = 'I型标尺，三序图')
points(db_rda, choices = 1:2, scaling = 1, display = 'sp', pch = 3, col = 'gray', cex = 1)
points(db_rda, choices = 1:2, scaling = 1, display = 'wa', pch = 19, col = c(rep('red', 9), rep('orange', 9), rep('green3', 9)), cex = 1)
text(db_rda, choices = 1:2, scaling = 1, display = 'cn', col = 'blue', cex = 0.8)

#############################
#使用原始物种多度数据执行 RDA，推荐先作个 Hellinger 转化
hel <- decostand(otu, method = 'hellinger')

#分别执行基于 Hellinger 转化后物种数据的 tb-RDA 以及与使用欧氏距离的 db-RDA
tb_rda_test <- rda(hel~., env, scale = FALSE)
db_rda_test <- capscale(hel~., env, distance = 'euclidean')

#可看到二者结果是一致的
par(mfrow = c(1, 2))
plot(tb_rda_test, scaling = 1, main = 'tb-RDA三序图')
plot(db_rda_test, scaling = 1, main = 'db-RDA三序图')

#但是基于卡方距离的 db-RDA，和 CCA 并不完全相同
cca_test <- cca(otu~., env, scale = FALSE)

chi <- decostand(otu, method = 'chi.square')	#先做卡方标准化，再计算欧氏距离，即可得卡方距离
db_rda_test <- capscale(chi~., env, distance = 'euclidean')

par(mfrow = c(2, 2))
plot(cca_test, scaling = 1, main = 'db-RDA三序图')
plot(db_rda_test, scaling = 1, main = 'CCA三序图')
plot(cca_test, display = c('wa', 'cn'), scaling = 1, main = 'db-RDA双序图')
plot(db_rda_test, display = c('wa', 'cn'), scaling = 1, main = 'CCA双序图')

#############################
##db-RDA 结果解读，以 capscale() 函数结果为例，简介
db_rda <- capscale(otu~., env, distance = 'bray', add = TRUE)

#查看统计结果信息，以 I 型标尺为例
db_rda.scaling1 <- summary(db_rda, scaling = 1)
db_rda.scaling1

#作图查看排序结果，详情 ?plot.cca
#三序图，包含 I 型标尺和 II 型标尺，样方坐标展示为使用物种加权计算的样方得分
par(mfrow = c(1, 2))

plot(db_rda, scaling = 1, main = 'I 型标尺', display = c('wa', 'sp', 'cn'))
rda_sp.scaling1 <- scores(db_rda, choices = 1:2, scaling = 1, display = 'sp')
#arrows(0, 0, rda_sp.scaling1[ ,1], rda_sp.scaling1[ ,2], length =  0, lty = 1, col = 'red')

plot(db_rda, scaling = 2, main = 'II 型标尺', display = c('wa', 'sp', 'cn'))
rda_sp.scaling2 <- scores(db_rda, choices = 1:2, scaling = 2, display = 'sp')
#arrows(0, 0, rda_sp.scaling2[ ,1], rda_sp.scaling2[ ,2], length =  0, lty = 1, col = 'red')

#隐藏物种，以 I 型标尺为例展示双序图
#比较分别使用物种加权计算的样方坐标以及拟合的样方坐标的差异
par(mfrow = c(1, 2))
plot(db_rda, scaling = 1, main = 'I 型标尺，加权', display = c('wa', 'cn'))
plot(db_rda, scaling = 1, main = 'I 型标尺，拟合', display = c('lc', 'cn'))

##RDA 结果提取
#scores() 提取排序得分（坐标），以 I 型标尺为例，前四轴为例
#使用物种加权和计算的样方得分
db_rda_site.scaling1 <- scores(db_rda, choices = 1:4, scaling = 1, display = 'wa')	
#物种变量（响应变量）得分
db_rda_sp.scaling1 <- scores(db_rda, choices = 1:4, scaling = 1, display = 'sp')
#环境变量（解释变量）得分
db_rda_env.scaling1 <- scores(db_rda, choices = 1:4, scaling = 1, display = 'bp')

#或者在 summary() 后提取，以 I 型标尺为例，前四轴为例
db_rda.scaling1 <- summary(db_rda, scaling = 1)
#使用物种加权和计算的样方得分
db_rda_site.scaling1 <- db_rda.scaling1$site[ ,1:4]
#物种
db_rda_sp.scaling1 <- db_rda.scaling1$species[ ,1:4]
#环境
db_rda_env.scaling1 <- db_rda.scaling1$biplot[ ,1:4]

#若需要输出在本地
#样方
write.table(data.frame(db_rda_site.scaling1), 'db_rda_site.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)
#物种
write.table(data.frame(db_rda_sp.scaling1), 'db_rda_sp.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)
#环境
write.table(data.frame(db_rda_env.scaling1), 'db_rda_env.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)

#不建议直接在原始数据集中提取，因为这样提取的坐标数值未经标尺缩放处理，不利于反映生物学问题
#db_rda$CCA$u[ ,1:4]
#db_rda$CCA$v[ ,1:4]
#db_rda$CCA$biplot[ ,1:4]

#RDA 基于多元回归，若对典范特征系数（即每个解释变量与每个约束轴之间的回归系数）感兴趣
#coef() 提取典范系数
rda_coef <- coef(db_rda)

#############################
##R2 校正
#RsquareAdj() 提取 R2，详情 ?RsquareAdj() 
r2 <- RsquareAdj(db_rda)
db_rda_noadj <- r2$r.squared	#原始 R2
db_rda_adj <- r2$adj.r.squared	#校正后的 R2

#关于约束轴承载的特征值或解释率，应当在 R2 校正后重新计算
db_rda_exp_adj <- db_rda_adj * db_rda$CCA$eig/sum(db_rda$CCA$eig)
db_rda_eig_adj <- db_rda_exp_adj * db_rda$tot.chi

##置换检验
#所有约束轴的置换检验，即全局检验，基于 999 次置换，详情 ?anova.cca
db_rda_test <- anova.cca(db_rda, permutations = 999)

#各约束轴逐一检验，基于 999 次置换
db_rda_test_axis <- anova.cca(db_rda, by = 'axis', permutations = 999)

#p 值校正（Bonferroni 为例）
db_rda_test_axis$`Pr(>F)` <- p.adjust(db_rda_test_axis$`Pr(>F)`, method = 'bonferroni')

##断棍模型和 Kaiser-Guttman 准则帮助确定重要的残差轴
pcoa_eig <- db_rda$CA$eig
n <- length(pcoa_eig)
bsm <- data.frame(j=seq(1:n), p = 0)
bsm$p[1] <- 1/n
for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
bsm$p <- 100*bsm$p/n
 
barplot(t(cbind(100 * pcoa_eig/sum(pcoa_eig), bsm$p[n:1])), beside = TRUE, main = '% 变差', col = c('orange', 'bisque'), las = 2)
abline(h = mean(100 * pcoa_eig/sum(pcoa_eig)), col = 'red')
legend('topright', c('% 特征根', '断棍模型', '平均特征根'), pch = c(15, 15, NA), col = c('orange', 'bisque', 'red'), lwd = c(NA, NA, 1), bty = 'n')

##变量选择
#计算方差膨胀因子，详情 ?vif.cca
vif.cca(db_rda)

#前向选择，以 ordiR2step() 的方法为例，基于 999 次置换检验，详情 ?ordiR2step
db_rda_forward_pr <- ordiR2step(capscale(otu~1, env, distance = 'bray', add = TRUE), scope = formula(db_rda), R2scope = TRUE, direction = 'forward', permutations = 999)

#以 db_rda 和 db_rda_forward_pr 为例，简要绘制双序图比较变量选择前后结果
par(mfrow = c(1, 2))
plot(db_rda, scaling = 1, main = '原始模型，I 型标尺', display = c('wa', 'cn'))
plot(db_rda_forward_pr, scaling = 1, main = '前向选择后，I 型标尺', display = c('wa', 'cn'))

#细节部分查看
summary(db_rda_forward_pr, scaling = 1)

#比较选择前后校正后 R2 的差异，详情 ?RsquareAdj
#可以看到变量选择后，尽管去除了很多环境变量，但总 R2 并未损失很多
RsquareAdj(db_rda)$adj.r.squared
RsquareAdj(db_rda_forward_pr)$adj.r.squared

#所有约束轴的全局检验，999 次置换，详情 ?anova.cca
db_rda_forward_pr_test <- anova.cca(db_rda_forward_pr, permutations = 999)

#各约束轴逐一检验，999 次置换
db_rda_forward_pr_test_axis <- anova.cca(db_rda_forward_pr, by = 'axis', permutations = 999)

#p 值校正（Bonferroni 为例）
db_rda_forward_pr_test_axis$`Pr(>F)` <- p.adjust(db_rda_forward_pr_test_axis$`Pr(>F)`, method = 'bonferroni')

##
#提取或输出变量选择后的排序坐标，以 I 型标尺为例，前两轴为例
#提取方式可参考上文，这里通过 scores() 提取

#使用物种加权和计算的样方得分
db_rda_forward_pr_site.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'wa')
write.table(data.frame(db_rda_forward_pr_site.scaling1), 'db_rda_forward_pr_site.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)
#物种变量（响应变量）得分
db_rda_forward_pr_sp.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'sp')
write.table(data.frame(db_rda_forward_pr_sp.scaling1), 'db_rda_forward_pr_sp.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)
#环境变量（解释变量）得分
db_rda_forward_pr_env.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'bp')
write.table(data.frame(db_rda_forward_pr_env.scaling1), 'db_rda_forward_pr_env.scaling1.txt', sep = '\t', col.names = NA, quote = FALSE)

##变差分解 varpart()
#以两组环境变量（DOC、AP+AK）为例，运行变差分解，详情 ?varpart
#输入 varpart() 的为已经计算好的 Bray-curtis 距离测度；参数 add=T 意为校正 PCoA 的负特征根（和上述 capscale(add=T) 对应）
db_rda_vp <- varpart(dis_bray, env['DOC'], env[c('AP', 'AK')], scale = FALSE, add = TRUE)
db_rda_vp

plot(db_rda_vp, digits = 2, Xnames = c('DOC', 'AP+AK'), bg = c('blue', 'red'))

#解释变差的置换检验，以 DOC 所能解释的全部变差为例；999 次置换
anova.cca(capscale(dis_bray~DOC, env, add = TRUE), permutations = 999)

#若考虑 DOC 单独解释的变差部分，需将其它变量作为协变量；999 次置换
anova.cca(capscale(dis_bray~DOC+Condition(AP+AK), env, add = TRUE), permutations = 999)

#############################
#以前向选择后的简约模型 db_rda_forward_pr 为例作图展示前两轴
#plot(db_rda_forward_pr) 方法可参考上文
#下面是 ggplot2

#提取样方和环境因子排序坐标，前两轴，I 型标尺
db_rda_forward_pr.scaling1 <- summary(db_rda_forward_pr, scaling = 1)
db_rda_forward_pr.site <- data.frame(db_rda_forward_pr.scaling1$sites)[1:2]
db_rda_forward_pr.env <- data.frame(db_rda_forward_pr.scaling1$biplot)[1:2]

#手动添加分组
db_rda_forward_pr.env$name <- rownames(db_rda_forward_pr.env)
db_rda_forward_pr.site$name <- rownames(db_rda_forward_pr.site)
db_rda_forward_pr.site$group <- c(rep('A', 9), rep('B', 9), rep('C', 9))

#计算校正 R2 后的约束轴解释率
exp_adj <- RsquareAdj(db_rda_forward_pr)$adj.r.squared * db_rda_forward_pr$CCA$eig/sum(db_rda_forward_pr$CCA$eig)
rda1_exp <- paste('RDA1:', round(exp_adj[1]*100, 2), '%')
rda2_exp <- paste('RDA1:', round(exp_adj[2]*100, 2), '%')

#ggplot2 作图
library(ggplot2)

p <- ggplot(db_rda_forward_pr.site, aes(CAP1, CAP2)) +
geom_point(aes(color = group)) +
scale_color_manual(values = c('red', 'orange', 'green3')) +
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) + 
labs(x = rda1_exp, y = rda2_exp) +
geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
geom_segment(data = db_rda_forward_pr.env, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'blue') +
geom_text(data = db_rda_forward_pr.env, aes(CAP1 * 1.1, CAP2 * 1.1, label = name), color = 'blue', size = 3)

p

