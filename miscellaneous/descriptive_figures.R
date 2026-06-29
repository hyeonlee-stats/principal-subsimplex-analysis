source('PSA_setup.R')
invisible(lapply(list.files('utils', pattern = '.R', full.names = T), source))

## ==== PCA failure ====
set.seed(1)
X = matrix(c(c(0.1,0.4,0.5) + rnorm(3*4, 0, 0.05),
             c(0.55,0.05,0.4) + rnorm(3*4, 0, 0.05)),
           byrow = T, ncol = 3)
X = matrix(c(c(0.1,0.45,0.45) + rnorm(3*4, 0, 0.05),
             c(0.5,0.05,0.45) + rnorm(3*4, 0, 0.05)),
           byrow = T, ncol = 3)
X = to_simplex(X)
X = rbind(X, matrix(c(0.01, 0.01, 0.98,
                      0.9, 0.05, 0.05),
                    byrow = T, ncol = 3))

X = as.data.frame(X)
colnames(X) = c('V1','V2','V3')

X.pca = princomp(X)
X.pc1 = rbind(X.pca$center - 2.5*X.pca$sdev[1]*X.pca$loadings[,1],
              X.pca$center + 2.5*X.pca$sdev[1]*X.pca$loadings[,1])
X.pc2 = rbind(X.pca$center - 2.5*X.pca$sdev[2]*X.pca$loadings[,2],
              X.pca$center + 2.5*X.pca$sdev[2]*X.pca$loadings[,2])
X.pc1 = to_2d(X.pc1)
X.pc2 = to_2d(X.pc2)

g = empty_tern(c('V1','V2','V3')) +
  geom_path(data = X.pc1, aes(x=V1, y=V2),
            arrow = arrow(length = unit(0.5, "cm")), col = 'forestgreen', linewidth = 1) +
  geom_path(data = X.pc2, aes(x=V1, y=V2),
            arrow = arrow(length = unit(0.5, "cm")), col = 'forestgreen', linewidth = 1) +
  geom_point(data = to_2d(matrix(X.pca$center, nrow=1)), col = 'black', size = 3) +
  geom_path(data = to_2d(rbind(X[9,], X.pca$center+X.pca$scores[9,1]*X.pca$loadings[,1])),
            col = 'black', linetype = 'dashed', linewidth = 1) +
  geom_path(data = to_2d(rbind(X[10,], X.pca$center+X.pca$scores[10,1]*X.pca$loadings[,1])),
            col = 'black', linetype = 'dashed', linewidth = 1) +
  geom_point(data = to_2d(X)[1:8,], col = 'blue', size = 3) +
  geom_point(data = to_2d(X)[9:10,], col = 'blue', size = 3) +
  geom_point(data = to_2d(matrix(rep(X.pca$center,2), byrow = T, nrow = 2) +
                          outer(X.pca$scores[9:10,1], X.pca$loadings[,1])), col = 'red', size = 3)
g$layers[[2]]$aes_params$size = 8

ggsave('auxiliary/Figures/PCA_failure.jpeg', g, width = 5.5, height = 5)

## ==== log-ratio example ====
### ==== data in original space ====
set.seed(1)
X = matrix(c(rep(c(0.15, 0.15, 0.7), 5),
              rep(c(0.3,0.1,0.6), 5)) + runif(3*10, 0, 0.2), ncol = 3, byrow = T)
X = to_simplex(X)
X = rbind(c(0.9, 0.05, 0.05),
           c(0.05, 0.05, 0.9),
           X)
X = as.data.frame(X)
colnames(X) = c('V1','V2','V3')

g1 = ggtern(X, aes(x=V1,y=V3,z=V2)) +
  geom_point(col = c('red','deepskyblue',rep('black',10))) +
  theme_bw()

### ==== modes of variation in clr space ====
X.clr = clr(X) %*% matrix(c(1,-1,0, -0.5,-0.5,1), byrow = F, ncol = 2)

pca.log = princomp(X.clr)
pc1.log = rbind(pca.log$center - 2.5*pca.log$sdev[1]*pca.log$loadings[,1],
              pca.log$center + 2*pca.log$sdev[1]*pca.log$loadings[,1])
pc2.log = rbind(pca.log$center - 7*pca.log$sdev[2]*pca.log$loadings[,2],
              pca.log$center + 7*pca.log$sdev[2]*pca.log$loadings[,2])

X.clr = as.data.frame(X.clr)
pc1.log = as.data.frame(pc1.log)
pc2.log = as.data.frame(pc2.log)
colnames(X.clr) = c('Z1','Z2')
colnames(pc1.log) = c('Z1','Z2')
colnames(pc2.log) = c('Z1','Z2')

g2 = ggplot(X.clr, aes(x=Z1,y=Z2)) +
  geom_point(col = c('red','deepskyblue',rep('black',10))) +
  geom_path(data = pc1.log, col = 'magenta') +
  geom_path(data = pc2.log, col = 'forestgreen') +
  theme_bw() +
  scale_x_continuous(limits = c(-1.5,3)) +
  scale_y_continuous(limits = c(-1.5,3))


### ==== modes of variation in original space ====
pca.clr = princomp(clr(X))
pc1.clr = as.data.frame(clrInv(outer(rep(1,2001),pca.clr$center) + outer((-1000:1000)/100,pca.clr$loadings[,1])))
pc2.clr = as.data.frame(clrInv(outer(rep(1,2001),pca.clr$center) + outer((-1000:1000)/100,pca.clr$loadings[,2])))
colnames(pc1.clr) = c('V1','V2','V3')
colnames(pc2.clr) = c('V1','V2','V3')

g3 = ggtern(X, aes(x=V1,y=V3,z=V2)) +
  geom_point(col = c('red','deepskyblue',rep('black',10))) +
  geom_path(data = pc1.clr, col = 'magenta') +
  geom_path(data = pc2.clr, col = 'forestgreen') +
  theme_bw()

g1 = ggtern::ggplot_gtable(ggtern::ggplot_build(g1))
g3 = ggtern::ggplot_gtable(ggtern::ggplot_build(g3))
plot_grid(g1, g2, g3, nrow = 1)
ggsave('auxiliary/Figures/logratio_pca.jpeg', width = 12, height = 4)

