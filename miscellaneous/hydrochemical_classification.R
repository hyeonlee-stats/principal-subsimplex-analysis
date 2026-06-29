source('PSA_setup.R')
library(e1071)

## read data
data(Hydrochem, package = 'compositions')
hc.label = Hydrochem$River
levels(hc.label) = list('Anoia' = 'Anoia', 'Cardener' = 'Cardner',
                        'Lower Llobregat' = 'LowerLLobregat', 'Upper Llobregat' = 'UpperLLobregat')
hc.X = to_simplex(as.matrix(Hydrochem[,6:19]))

## get scores
hc.pca = comp_pca(hc.X)
hc.power_pca = comp_power_pca(hc.X, 0.5)
hc.apca = comp_apca(hc.X)
hc.psas = psa(hc.X, 's')
hc.psao = psa(hc.X, 'o')

## comparison
set.seed(1)
L = 10

hc.pca.accuracy = rep(NA, L)
hc.pca.df = data.frame(label = hc.label) %>%
  cbind(as.data.frame(hc.pca$scores))
for(i in 1:L){
  svm.fit = tune.svm(label ~., data = hc.pca.df[,1:(i+1)],
                     kernel = 'radial',
                     gamma = 2**seq(-4,4,length.out=17),
                     cost = 1,
                     tunecontrol = tune.control(sampling = "cross", cross = 5))

  hc.pca.accuracy[i] = 1-svm.fit$best.performance
}

hc.power_pca.accuracy = rep(NA, L)
hc.power_pca.df = data.frame(label = hc.label) %>%
  cbind(as.data.frame(hc.power_pca$scores))
for(i in 1:L){
  svm.fit = tune.svm(label ~., data = hc.power_pca.df[,1:(i+1)],
                     kernel = 'radial',
                     gamma = 2**seq(-4,4,length.out=17),
                     cost = 1,
                     tunecontrol = tune.control(sampling = "cross", cross = 5))

  hc.power_pca.accuracy[i] = 1-svm.fit$best.performance
}

hc.apca.accuracy = rep(NA, L)
hc.apca.df = data.frame(label = hc.label) %>%
  cbind(as.data.frame(hc.apca$scores))
for(i in 1:L){
  svm.fit = tune.svm(label ~., data = hc.apca.df[,1:(i+1)],
                     kernel = 'radial',
                     gamma = 2**seq(-4,4,length.out=17),
                     cost = 1,
                     tunecontrol = tune.control(sampling = "cross", cross = 5))

  hc.apca.accuracy[i] = 1-svm.fit$best.performance
}

hc.psas.accuracy = rep(NA, L)
for(i in 1:L){
  df = as.data.frame(hc.psas$Xhat_reduced[[i+1]]) %>%
    mutate(label = hc.label)
  svm.fit = tune.svm(label ~., data = df,
                     kernel = 'radial',
                     gamma = 2**seq(-4,4,length.out=17),
                     cost = 1,
                     tunecontrol = tune.control(sampling = "cross", cross = 5))

  hc.psas.accuracy[i] = 1-svm.fit$best.performance
}

hc.psao.accuracy = rep(NA, L)
for(i in 1:L){
  df = as.data.frame(hc.psao$Xhat_reduced[[i+1]]) %>%
    mutate(label = hc.label)
  svm.fit = tune.svm(label ~., data = df,
                     kernel = 'radial',
                     gamma = 2**seq(-4,4,length.out=17),
                     cost = 1,
                     tunecontrol = tune.control(sampling = "cross", cross = 5))

  hc.psao.accuracy[i] = 1-svm.fit$best.performance
}


################################################################################

accuracy = as.data.frame(rbind(hc.pca.accuracy,
                                   hc.power_pca.accuracy,
                                   hc.psas.accuracy,
                                   hc.psao.accuracy,
                                   hc.apca.accuracy)) %>%
  mutate(Method = factor(c('PCA','Power transform','PSAS','PSAO','Log-ratio'),
                         levels = c('PSAS','PSAO','PCA','Power transform','Log-ratio')))

accuracy %>%
  melt(id.vars = 'Method') %>%
  mutate(variable = as.numeric(as.factor(variable))) %>%
  ggplot(aes(x=variable, y=value, group = Method, color = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  scale_x_continuous(name = 'Rank', breaks = 1:L) +
  scale_y_continuous(name = 'Accuracy', limits = c(0.4,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  scale_color_manual(values = c('red','orange','blue','dodgerblue','darkgray'))
ggsave('auxiliary/Figures/hydrochemical_accuracy.jpeg', width = 8, height = 4)


parallel_coord(hc.X, hc.label) +
  theme(legend.position='right') +
  scale_color_discrete(name = 'Tributary')
ggsave('auxiliary/Figures/hydrochemical_parallel_coord.jpeg', width = 8, height = 4)
