source('PSA_setup.R')

set.seed(1)
X = rbind(make_cluster(c(0.1, 0.8, 0.1), 10, 0.1^2),
          make_cluster(c(0.8, 0.1, 0.1), 10, 0.1^2),
          make_cluster(c(0.45, 0.45, 0.1), 30, 0.1^2))
df = as.data.frame(X) %>%
  mutate(label = factor(ifelse(V1>0.45, '1', '2')))
X.pca = comp_pca(X)
X.power_pca = comp_power_pca(X, 0.5)
X.psas = psa(X, 's')
X.psao = psa(X, 'o')
X.apca = comp_apca(X)

get_error_rate <- function(y_true){
  n = length(y_true)
  error_rate = rep(0., n)
  for(i in 1:n){
    y1 = y_true[1:i]
    y2 = y_true[(i+1):n]
    e = (sum(y1!='1') + sum(y2!='2'))/n
    error_rate[i] = min(e, 1-e)
  }
  return(error_rate)
}

c(data.frame(PC1 = X.psas$scores[,1],
             label = df$label) %>%
    arrange(PC1) %>%
    pull(label) %>%
    get_error_rate() %>%
    min(na.rm = T),
  data.frame(PC1 = X.psao$scores[,1],
             label = df$label) %>%
    arrange(PC1) %>%
    pull(label) %>%
    get_error_rate() %>%
    min(na.rm = T),
  data.frame(PC1 = X.pca$scores[,1],
             label = df$label) %>%
    arrange(PC1) %>%
    pull(label) %>%
    get_error_rate() %>%
    min(na.rm = T),
  data.frame(PC1 = X.power_pca$scores[,1],
             label = df$label) %>%
    arrange(PC1) %>%
    pull(label) %>%
    get_error_rate() %>%
    min(na.rm = T),
  data.frame(PC1 = X.apca$scores[,1],
             label = df$label) %>%
    arrange(PC1) %>%
    pull(label) %>%
    get_error_rate() %>%
    min(na.rm = T))

plot_grid(ternary_pc(X, df$label, type = 'data'),
          ternary_pc(X, df$label, type = 'psas', X.psas) +
            annotate('text',x=0.7,y=1.5,label='error rate: 0.00',size=6),
          ternary_pc(X, df$label, type = 'psao', X.psao) +
            annotate('text',x=0.7,y=1.5,label='error rate: 0.08',size=6),
          ternary_pc(X, df$label, type = 'pca', X.pca) +
            annotate('text',x=0.7,y=1.5,label='error rate: 0.08',size=6),
          ternary_pc(X, df$label, type = 'power', X.power_pca) +
            annotate('text',x=0.7,y=1.5,label='error rate: 0.06',size=6),
          ternary_pc(X, df$label, type = 'logratio', X.apca) +
            annotate('text',x=0.7,y=1.5,label='error rate: 0.20',size=6),
          nrow = 2)
ggsave('auxiliary/Figures/classification_example.jpeg', width = 12, height = 8)
