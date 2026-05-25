source('PSA_setup.R')

data.path = 'Diatom/Data/'
image.path = 'Diatom/Figures/'
dir.create(data.path, showWarnings = T)
dir.create(image.path, showWarnings = T)

## ==== diatom preparation ====
### ==== read data ====
diatom.df = read.csv(paste0(data.path, 'diatom.csv'))
diatom.info = read.csv(paste0(data.path,'diatom codes.csv'))

diatom.X = as.matrix(diatom.df[,-1])
colnames(diatom.X) = diatom.info$code
diatom.X = diatom.X[,names(sort(colMeans(diatom.X), decreasing = T))] # sort by mean proportion

rownames(diatom.info) = diatom.info$code
diatom.info$class[diatom.info$class == ''] = 'others'
diatom.info$class = factor(diatom.info$class, levels = c('warm water','open ocean','sea ice','others'))
diatom.info$color = c('hotpink','green3','dodgerblue','black')[diatom.info$class]

### ==== get scores ====
diatom.res = list(X = diatom.X, info = diatom.info)
# system.time({diatom.res$psas = psa(diatom.X, 's')})
# saveRDS(diatom.res$psas, paste0(data.path, 'diatom_psas.rds'))
# system.time({diatom.res$psao = psa(diatom.X, 'o')})
# saveRDS(diatom.res$psao, paste0(data.path, 'diatom_psao.rds'))
diatom.res$psas = readRDS(paste0(data.path, 'diatom_psas.rds'))
diatom.res$psao = readRDS(paste0(data.path, 'diatom_psao.rds'))
diatom.res$pca = comp_pca(diatom.X)
diatom.res$power_pca = comp_power_pca(diatom.X, 0.5)
diatom.res$apca = comp_apca(diatom.X)

# flip direction of loadings and scores in accordance with PSA
diatom.res$pca = flip_loading(diatom.res$pca, c(1,2,3))
diatom.res$power_pca = flip_loading(diatom.res$power_pca, c(1,2,4))
diatom.res$apca = flip_loading(diatom.res$apca, c(1,2))

## ==== descriptive figures ====
### ==== marker setup ====
diatom.colors = list(colors = c('blue','magenta','red','orange','green'),
                     values = rescale(c(56.53,65.04,67.03,69.74,83.99), to = c(0,1)))
g = ggplot(diatom.df, aes(x=Depth, y=Depth, col=Depth)) +
  scale_color_gradientn(colors = diatom.colors$colors,
                        values = diatom.colors$values,
                        guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  geom_point() +
  theme(legend.key.size = unit(14, 'points'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.spacing.y = unit(1,'points'))
diatom.legend = get_legend(g)
diatom.col = ggplot_build(g)$data[[1]]$colour
diatom.shape = c(rep(1,35),rep(16,4),rep(1,32))
diatom.size = c(rep(1,35),rep(2,4),rep(1,32))

### ==== heatmap ====
diatom.X.log = log10(replace_zero(diatom.X))

heatmap.df = as.data.frame(diatom.X.log) %>%
  cbind(data.frame(Depth = diatom.df$Depth)) %>%
  mutate(Depth_idx = rank(-Depth)) %>%
  select(-Depth) %>%
  melt(id.vars = 'Depth_idx') %>%
  mutate(variable = factor(variable, levels = colnames(diatom.X)))

g2 = ggplot(heatmap.df, aes(x = variable, y = Depth_idx, fill = value)) +
  theme_minimal() +
  geom_tile() +
  scale_fill_gradientn(colors = c('white','black','black'),
                       values = rescale(c(max(diatom.X.log), max(apply(diatom.X.log, 2, min)), min(diatom.X.log))),
                       name = 'log10(Proportion)') +
  labs(y = 'Depth (Index)') +
  scale_y_continuous(expand = c(0,0,0,0)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_discrete(labels = diatom.info[as.character(colnames(diatom.X)),'code']) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'],
                                   angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_text(hjust=0.5))

g1 = data.frame(variable = 'Depth',
                Depth = diatom.df$Depth) %>%
  ggplot(aes(x = variable, y = rank(Depth), fill = Depth)) +
  theme_minimal() +
  geom_tile() +
  scale_fill_gradientn(colors = diatom.colors$colors,
                       values = diatom.colors$values,
                       guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  scale_y_continuous(trans = 'reverse',
                     labels = diatom.df$Depth[c(1,18,36,54,71)], breaks = c(1,18,36,54,71),
                     expand = c(0,0,0,0)) +
  labs(x = NULL, y = 'Depth') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.position = 'none')

g = plot_grid(g1, g2, nrow = 1, rel_widths = c(1,10), align = 'h', axis = 'tb')
ggsave('1. diatom_heatmap.jpeg', g, path = image.path, width = 12, height = 4)


### ==== mean proportion versus log scale variance ====
g = data.frame(mean_prop = colMeans(diatom.X),
               sd_log10 = apply(diatom.X.log, 2, sd),
               label = diatom.info$class) %>%
  ggplot(aes(x = mean_prop, y = sd_log10, col = label, shape = label)) +
  theme_bw() +
  geom_point() +
  scale_color_manual(values = c('hotpink','green3','dodgerblue','black'),
                     name = 'Type') +
  scale_shape_manual(values = c(16,16,16,1),
                     name = 'Type') +
  labs(x = 'Mean proportion', y = 'Standard deviation (log10)')
ggsave('1. diatom_mean_proportion_versus_variance.jpeg', g, path = image.path, width = 6, height = 3)

### ==== parallel coordinate plot ====
g = parallel_coord(diatom.X, diatom.df$Depth, diatom.info[as.character(colnames(diatom.X)),'code']) +
  scale_color_gradientn(colors = diatom.colors$colors,
                        values = diatom.colors$values,
                        name = 'Depth',
                        guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  labs(x = '', y = 'Proportion') +
  theme(legend.position = 'right', legend.key.size = unit(10, 'points'), legend.title = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.key.size = unit(12, 'points'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing.y = unit(1,'points')) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color']))

ggsave('1. Parallel_plot_raw.jpeg', g, path = image.path, width = 12, height = 6)

### ==== depth distribution ====
n = length(diatom.df$Depth)
diatom.depth.df = data.frame(Depth = diatom.df$Depth,
                             y = rev((1:n)/n-1/2/n))
g = ggplot(diatom.depth.df, aes(x = Depth, y = y, col = Depth)) +
  scale_x_continuous(trans = 'reverse') +
  geom_point(shape = 1.5, size = 2, stroke = 1) +
  labs(y='') +
  theme_bw() +
  scale_color_gradientn(colors = diatom.colors$colors,
                        values = diatom.colors$values,
                        name = 'Depth',
                        guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  theme(legend.key.size = unit(8, 'points'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.spacing.y = unit(1,'points'),
        plot.margin = unit(rep(15,4),'points'),
        axis.title.x = element_text(size = 8))

ggsave('1. Depth_distribution.jpeg', g, path = image.path, width = 8, height = 3)

### ==== percent nonzero ====
df = data.frame(variable = colnames(diatom.X),
                percent.nonzero = colSums(diatom.X > 0)/nrow(diatom.X))
df$variable = factor(df$variable, levels = df$variable)
g = ggplot(df, aes(x = variable, y = percent.nonzero)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = '', y = 'Proportion Nonzero') +
  geom_hline(yintercept = c(0.25,0.5,0.75), col = 'red') +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'])) +
  scale_x_discrete(labels = diatom.info[as.character(colnames(diatom.X)),'code'])

ggsave('1. Percent_nonzero.jpeg', g, path = image.path, width = 8, height = 4)

## ==== Comparison ====
### ==== create all panels ====
g.psas = line_gpairs(cbind(diatom.res$psas$scores[,1:4], data.frame(Depth = diatom.df$Depth)),
                     col.point = diatom.col, col.line = diatom.col,
                     shape = diatom.shape, size = diatom.size)
g.psao = line_gpairs(cbind(diatom.res$psao$scores[,1:4], data.frame(Depth = diatom.df$Depth)),
                     col.point = diatom.col, col.line = diatom.col,
                     shape = diatom.shape, size = diatom.size)
g.pca = line_gpairs(cbind(diatom.res$pca$scores[,1:4], data.frame(Depth = diatom.df$Depth)),
                    col.point = diatom.col, col.line = diatom.col,
                    shape = diatom.shape, size = diatom.size)
g.power_pca = line_gpairs(cbind(diatom.res$power_pca$scores[,1:4], data.frame(Depth = diatom.df$Depth)),
                          col.point = diatom.col, col.line = diatom.col,
                          shape = diatom.shape, size = diatom.size)
g.apca = line_gpairs(cbind(diatom.res$apca$scores[,1:4], data.frame(Depth = diatom.df$Depth)),
                     col.point = diatom.col, col.line = diatom.col,
                     shape = diatom.shape, size = diatom.size)

### ==== score 1 vs score 2 ====
g = plot_grid(diatom.legend,
              getPlot(g.psas, 2, 1) + ggtitle('PSA-S') +
                scale_x_continuous(limits=c(-0.4,0.63)) + scale_y_continuous(limits=c(-0.4,0.63)),
              getPlot(g.psao, 2, 1) + ggtitle('PSA-O') +
                scale_x_continuous(limits=c(-0.52,0.81)) + scale_y_continuous(limits=c(-0.52,0.81)),
              getPlot(g.pca, 2, 1) + ggtitle('PCA') +
                scale_x_continuous(limits=c(-0.21,0.23)) + scale_y_continuous(limits=c(-0.21,0.23)),
              getPlot(g.power_pca, 2, 1) + ggtitle('Power Transform PCA') +
                scale_x_continuous(limits=c(-0.12,0.2)) + scale_y_continuous(limits=c(-0.12,0.2)),
              getPlot(g.apca, 2, 1) + ggtitle('Log-ratio PCA') +
                scale_x_continuous(limits=c(-7.3,7.3)) + scale_y_continuous(limits=c(-7.3,7.3)),
              nrow = 2, align='hv', axis='tblr')
ggsave('2. Score_1vs2.jpeg', g, path = image.path, width = 10.5, height = 7)

### ==== PSA-O scores ====
g = plot_grid(getPlot(g.psao, 5, 1),
              getPlot(g.psao, 5, 2),
              getPlot(g.psao, 5, 3),
              getPlot(g.psao, 5, 4),
              diatom.legend,
              nrow = 1, rel_widths = c(2,2,2,2,1))
ggsave('2. Score_psao_diag.jpeg', g, path = image.path, width = 14.5, height = 3)

### ==== save matrices ====
add_gap <- function(g){
  gg = ggmatrix_gtable(gpairs_lower(g))
  gg$heights[[15]] = unit(0.4, 'null')
  gg
}

g.psas = add_gap(g.psas)
g.psao = add_gap(g.psao)
g.pca = add_gap(g.pca)
g.power_pca = add_gap(g.power_pca)
g.apca = add_gap(g.apca)

g = plot_grid(diatom.legend,
              g.psas, g.psao, g.pca, g.power_pca, g.apca, nrow = 2)
ggsave('2. Score_plot_matrix.jpeg', g, path = image.path, width = 18, height = 12)

g.psas = plot_grid(g.psas, diatom.legend, nrow = 1, rel_widths = c(5,1))
g.psao = plot_grid(g.psao, diatom.legend, nrow = 1, rel_widths = c(5,1))
g.pca = plot_grid(g.pca, diatom.legend, nrow = 1, rel_widths = c(5,1))
g.power_pca = plot_grid(g.power_pca, diatom.legend, nrow = 1, rel_widths = c(5,1))
g.apca = plot_grid(g.apca, diatom.legend, nrow = 1, rel_widths = c(5,1))

ggsave('2-1. Score_psas.jpeg', g.psas, path = image.path, width = 9*1.2, height = 9)
ggsave('2-2. Score_psao.jpeg', g.psao, path = image.path, width = 9*1.2, height = 9)
ggsave('2-3. Score_pca.jpeg', g.pca, path = image.path, width = 9*1.2, height = 9)
ggsave('2-4. Score_power_pca.jpeg', g.power_pca, path = image.path, width = 9*1.2, height = 9)
ggsave('2-5. Score_apca.jpeg', g.apca, path = image.path, width = 9*1.2, height = 9)

### ==== loading bar plots ====
plot_loadings_diatom <- function(V, k = 4, max.k = 12){
  ls = list()
  for(i in 1:k){
    gg = plot_vertex(V[,i], max.k) +
      ggtitle(colnames(V)[i]) +
      theme(axis.text.y = element_text(size=12))
    ls[[i]] = gg + theme(axis.text.y = element_text(color = diatom.info[as.character(gg$data$variable),'color']))
  }
  g = cowplot::plot_grid(plotlist = ls, nrow = 1, align = 'v', axis = 'l')
  return(g)
}

g = plot_loadings_diatom(diatom.res$psas$loadings)
ggsave('3. loading_bar_psas.jpeg', g, path = image.path, width = 16, height = 4)
g = plot_loadings_diatom(diatom.res$psao$loadings)
ggsave('3. loading_bar_psao.jpeg', g, path = image.path, width = 16, height = 4)

g = plot_grid(ggdraw() + draw_label('PSA-S'), plot_loadings_diatom(diatom.res$psas$loadings),
              ggdraw() + draw_label('PSA-O'), plot_loadings_diatom(diatom.res$psao$loadings),
              ggdraw() + draw_label('PCA'), plot_loadings_diatom(diatom.res$pca$loadings),
              ggdraw() + draw_label('Power Transform PCA'), plot_loadings_diatom(diatom.res$power_pca$loadings),
              ggdraw() + draw_label('Log-ratio PCA'), plot_loadings_diatom(diatom.res$apca$loadings),
          ncol=2, align='v', axis='lrtb', rel_widths = c(1,6))
ggsave('3. loading_bar_all.jpeg', g, path = image.path, width = 16, height = 16)

### ==== loading parallel coord ====
g.axis = loading_parallel(diatom.res$psas$loadings[,1],
                          diatom.info[as.character(colnames(diatom.X)),'code']) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color']))

g = plot_grid(plot_grid(ggdraw() + draw_label('Comp.1'), ggdraw() + draw_label('Comp.2'),
                        ggdraw() + draw_label('Comp.3'), ggdraw() + draw_label('Comp.4'),
                        ggdraw(),
                        ncol = 1, rel_heights = c(1,1,1,1,3)),
              plot_grid(loading_parallel(diatom.res$psas$loadings[,1]) +
                          theme(axis.text.x = element_blank()),
                        loading_parallel(diatom.res$psas$loadings[,2]) +
                          theme(axis.text.x = element_blank()),
                        loading_parallel(diatom.res$psas$loadings[,3]) +
                          theme(axis.text.x = element_blank()),
                        loading_parallel(diatom.res$psas$loadings[,4]) +
                          theme(axis.text.x = element_blank()),
                        ggdraw(get_x_axis(g.axis)),
                        ncol = 1, rel_heights = c(1,1,1,1,2), align = 'v', axis = 'lr'),
              nrow = 1, rel_widths = c(1,9))

ggsave('4. loading_parallel_psas.jpeg', g, path = image.path, width = 12, height = 7)

g = plot_grid(plot_grid(ggdraw() + draw_label('Comp.1'), ggdraw() + draw_label('Comp.2'),
                        ggdraw() + draw_label('Comp.3'), ggdraw() + draw_label('Comp.4'),
                        ggdraw(),
                        ncol = 1, rel_heights = c(1,1,1,1,3)),
              plot_grid(loading_parallel(diatom.res$psao$loadings[,1]) +
                          theme(axis.text.x = element_blank()),
                        loading_parallel(diatom.res$psao$loadings[,2]) +
                          theme(axis.text.x = element_blank()),
                        loading_parallel(diatom.res$psao$loadings[,3]) +
                          theme(axis.text.x = element_blank()),
                        loading_parallel(diatom.res$psao$loadings[,4]) +
                          theme(axis.text.x = element_blank()),
                        ggdraw(get_x_axis(g.axis)),
                        ncol = 1, rel_heights = c(1,1,1,1,2), align = 'v', axis = 'lr'),
              nrow = 1, rel_widths = c(1,9))

ggsave('4. loading_parallel_psao.jpeg', g, path = image.path, width = 12, height = 7)

## ==== ternary ====
### ==== ternary plot ====
psas.tern = as.data.frame(diatom.res$psas$Xhat_reduced$`r=2`) %>%
  mutate(Depth = diatom.df$Depth) %>%
  ggtern(aes(x=V1, y=V3, z=V2)) +
  geom_path(col = diatom.col, linewidth = 0.5) +
  geom_point(col = diatom.col, shape = diatom.shape, size = diatom.size*4/3, stroke = 1) +
  theme_bw() +
  theme(legend.position = 'none')

psas.tern = plot_grid(ggtern::ggplot_gtable(ggtern::ggplot_build(psas.tern)),
                      diatom.legend,
                      rel_widths = c(3,1))
ggsave('3-1. ternary_psas.jpeg', psas.tern, path = image.path, width = 8, height = 5)

psao.tern = as.data.frame(diatom.res$psao$Xhat_reduced$`r=2`) %>%
  mutate(Depth = diatom.df$Depth) %>%
  ggtern(aes(x=V1, y=V3, z=V2)) +
  geom_path(col = diatom.col, linewidth = 0.5) +
  geom_point(col = diatom.col, shape = diatom.shape, size = diatom.size*4/3, stroke = 1) +
  theme_bw() +
  theme(legend.position = 'none')

psao.tern = plot_grid(ggtern::ggplot_gtable(ggtern::ggplot_build(psao.tern)),
                      diatom.legend,
                      rel_widths = c(3,1))
ggsave('3-2. ternary_psao.jpeg', psao.tern, path = image.path, width = 8, height = 5)

### ==== ternary vertices ====
g = plot_loadings_diatom(diatom.res$psas$Vhat$`r=2`, 3)
ggsave('3-3. ternary_vertex_psas.jpeg', g, path = image.path, width = 12, height = 4)
g = plot_loadings_diatom(diatom.res$psao$Vhat$`r=2`, 3)
ggsave('3-4. ternary_vertex_psao.jpeg', g, path = image.path, width = 12, height = 4)


## ==== Variance explained & RSS ====
### ==== variance explained ====
var_exp = rbind(c(diatom.res$psas$RSS,0),
            c(diatom.res$psao$RSS,0),
            diatom.res$pca$RSS,
            diatom.res$power_pca$RSS,
            diatom.res$apca$RSS)
rownames(var_exp) = c('PSA-S','PSA-O','PCA','Power transform','Log-ratio')
g = plot_variance_explained(var_exp, 6)
ggsave('variance_explained.jpeg', g, path = image.path, width = 8, height = 4)

### ==== RSS ====
rss.mat = matrix(NA, 6, 6)
rownames(rss.mat) = c('PSA-S','PSA-O','PCA','Power transform','Log-ratio','PCA (projected)')
colnames(rss.mat) = paste0('PC',1:6)

for(r in 1:6){
  r.str = paste0('r=',r)
  rss.mat[,r] = c(sum(rowSums((diatom.res$psas$Xhat[[r.str]]-diatom.res$X)^2)),
                  sum(rowSums((diatom.res$psao$Xhat[[r.str]]-diatom.res$X)^2)),
                  sum(rowSums((diatom.res$pca$Xhat[[r.str]]-diatom.res$X)^2)),
                  sum(rowSums((diatom.res$power_pca$Xhat[[r.str]]-diatom.res$X)^2)),
                  sum(rowSums((diatom.res$apca$Xhat[[r.str]]-diatom.res$X)^2)),
                  sum(rowSums((to_simplex(diatom.res$pca$Xhat[[r.str]])-diatom.res$X)^2)))
}

g = rss.mat %>%
  reshape2::melt() %>%
  setNames(c('Method','Rank','value')) %>%
  mutate(Method = factor(Method, levels = c('PSA-S','PSA-O','PCA','Power transform','Log-ratio','PCA (projected)'))) %>%
  ggplot(aes(x = Rank, y = value, group = Method, shape = Method, linetype = Method)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  scale_linetype_manual(values = c(rep('solid',5), 'dashed')) +
  scale_shape_manual(values = c(16,17,15,3,7,15)) +
  ylab('Residual sums of squares')
ggsave('diatom_rss.jpeg', g, path = image.path, width = 8, height = 4)


## ==== biplots ====
draw_biplot <- function(res){
  x = res$scores[,1:2]
  y = res$loadings[,1:2]
  y = y/sqrt(mean(rowSums(y**2))/mean(rowSums(x**2)))/3

  dfx = data.frame(x = x[,1], y = x[,2], Depth = diatom.df$Depth)
  dfy = data.frame(x = 0, y = 0, xend = y[,1], yend = y[,2],
                   size = sqrt(rowSums(y[,1:2]**2)),
                   label=diatom.info[rownames(y),]$short_code,
                   color=factor(diatom.info[rownames(y),]$color))
  dfy = dfy %>% arrange(desc(size)) %>% head(12)
  dfy2 = dfy
  yend = dfy2$yend
  dfy2$yend[abs(yend)<0.001] = (sum(abs(yend)<0.001):1)*(max(yend)-min(yend))/20

  ggplot() +
    theme_bw() +
    geom_point(data=dfx, aes(x=x,y=y,col=Depth)) +
    geom_path(data=dfx, aes(x=x,y=y,col=Depth)) +
    scale_color_gradientn(colors = c('blue','magenta','red','orange','green'),
                          values = rescale(c(56.53,65.04,67.03,69.74,83.99), to = c(0,1)),
                          name = 'Depth',
                          guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
    geom_segment(data=dfy, aes(x=x,y=y,xend=xend,yend=yend),
                 arrow = arrow(length=unit(5,'points')), linewidth=0.2, color = 'black') +
    new_scale_color() +
    geom_text(data=dfy2, aes(x=xend, y=yend, label=label, col=color), size=3) +
    scale_color_manual(values = levels(dfy$color)) +
    labs(x='Comp.1', y='Comp.2') +
    theme(legend.position = 'none',
          plot.title=element_text(hjust=0.5))
}

g = plot_grid(draw_biplot(diatom.res$pca) + ggtitle('PCA') +
                scale_x_continuous(limits=c(-0.21,0.23)) + scale_y_continuous(limits=c(-0.21,0.23)),
              draw_biplot(diatom.res$power_pca) + ggtitle('Power Transform PCA') +
                scale_x_continuous(limits=c(-0.12,0.2)) + scale_y_continuous(limits=c(-0.12,0.2)),
              draw_biplot(diatom.res$apca) + ggtitle('Log-ratio PCA') +
                scale_x_continuous(limits=c(-7.3,7.3)) + scale_y_continuous(limits=c(-7.3,7.3)),
              diatom.legend,
              draw_biplot(diatom.res$psas) + ggtitle('PSA-S') +
                scale_x_continuous(limits=c(-0.4,0.63)) + scale_y_continuous(limits=c(-0.4,0.63)),
              draw_biplot(diatom.res$psao) + ggtitle('PSA-O') +
                scale_x_continuous(limits=c(-0.52,0.81)) + scale_y_continuous(limits=c(-0.52,0.81)),
              nrow = 3, align='hv', axis='tblr', byrow = F)
ggsave('score_biplot.jpeg', g, path = image.path, width = 7, height = 10.5)


## ==== PCA before projection ====
g = parallel_coord(diatom.res$pca$Xhat$`r=1`, diatom.df$Depth, diatom.info[as.character(colnames(diatom.X)),'code']) +
  scale_color_gradientn(colors = diatom.colors$colors,
                        values = diatom.colors$values,
                        name = 'Depth',
                        guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  labs(x = '', y = 'Proportion') +
  theme(legend.position = 'right', legend.key.size = unit(10, 'points'), legend.title = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.key.size = unit(12, 'points'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing.y = unit(1,'points')) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'])) +
  geom_hline(yintercept = 0, col='black', linewidth=1) +
  scale_y_continuous(limits = c(-0.1,0.2))
ggsave('parallel_coord_pca_r1.jpeg', g, path = image.path, width = 12, height = 6)

g = parallel_coord(to_simplex(diatom.res$pca$Xhat$`r=1`), diatom.df$Depth, diatom.info[as.character(colnames(diatom.X)),'code']) +
  scale_color_gradientn(colors = diatom.colors$colors,
                        values = diatom.colors$values,
                        name = 'Depth',
                        guide = guide_colorbar(direction = "vertical", reverse = TRUE)) +
  labs(x = '', y = 'Proportion') +
  theme(legend.position = 'right', legend.key.size = unit(10, 'points'), legend.title = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.key.size = unit(12, 'points'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing.y = unit(1,'points')) +
  theme(axis.text.x = element_text(color = diatom.info[as.character(colnames(diatom.X)),'color'])) +
  geom_hline(yintercept = 0, col='black', linewidth=1) +
  scale_y_continuous(limits = c(-0.1,0.2))
ggsave('parallel_coord_pca_r1_projected.jpeg', g, path = image.path, width = 12, height = 6)

