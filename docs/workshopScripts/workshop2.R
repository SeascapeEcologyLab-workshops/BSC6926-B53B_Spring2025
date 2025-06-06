#' """ Workshop 2: Trophic niche dynamics
#'     @author: BSC6926-B53B
#'     date: 3/21/25"""

library(tidyverse)
library(hypervolume)

#Function to make random points of n length----
# data from random sample with mean and sd but 
# can be generated between a high and low value if end_points = T
# *** Note chose either column names or column numbers for ID_rows and names
# either work but must be the same
# df = dataframe or tibble with each row containing 
#        unique entry for making random points 
# ID_rows = vector of column names or numbers with id information
# names = column name or number of name of measure variables
# mean = column name or column number of df with mean 
# sd = column name or column number of df with sd 
# n = number of points to randomly generate
# z_score = T or F, if T z-scores values
# end_points = T or F for if random points need to be generated between
#        a high and low end point (e.g. 5% and 95% interval)
#        low and high required if end_points = T
# low = column name or column number of df with lower bound to sample in
# high = column name or column number of df with upper bound to sample in
HVvalues = function(df, ID_rows, names, mean, sd, n, z_score = F,
                    end_points = F, low = NULL, high = NULL){
  require(tidyverse)
  require(truncnorm)
  
  # check to see if information is needed to restrict where points are
  if (end_points){
    if (is_empty(df[,low]) | is_empty(df[,high])){
      return(cat('Warning: low and/or high columns not specified \n
                  Specific and run again or end_points = F \n'))
    }
  }
  
  # check to see if there are more 
  if(T %in% duplicated(df[,c(ID_rows,names)])){
    return(cat('Warning: some of the rows contain duplicated information \n
                make sure data is correct \n'))
  }
  
  # rename variables to make code work
  if (is.numeric(mean)){
    names(df)[mean] = 'Mean'
  }else {
    df = df |>  rename(Mean = mean)
  }
  
  if (is.numeric(sd)){
    names(df)[sd] = 'SD'
  }else {
    df = df |>  rename(SD = sd)
  }
  
  if (end_points){
    if (is.numeric(low)){
      names(df)[low] = 'lower'
    }else {
      df = df |>  rename(lower = low)
    }
    
    if (is.numeric(high)){
      names(df)[high] = 'upper'
    }else {
      df = df |>  rename(upper = high)
    }
  }
  
  # make sure the names is not numeric 
  if (is.numeric(names)){
    names = names(df)[names]
  }
  
  labs = unique(as.vector(df[,names])[[1]])
  # generate random points within bounds
  if(end_points){
    
    df_tot = df |> slice(rep(1:n(), each=n))|> 
      mutate(point = 
               truncnorm::rtruncnorm(1, a = lower, b = upper,
                                     mean = Mean, sd = SD)) |> 
      ungroup() |> 
      mutate(num = rep(1:n, times=nrow(df))) |>
      dplyr::select(-Mean, -SD, -lower, -upper)|>
      pivot_wider(names_from = names, values_from = point)|> 
      dplyr::select(-num)
  }else {
    # generate random points outside of bounds
    df_tot = df |> slice(rep(1:n(), each=n))|>
      mutate(point = 
               truncnorm::rtruncnorm(1, mean = Mean, sd = SD)) |> 
      ungroup() |> 
      mutate(num = rep(1:n, times=nrow(df))) |>
      dplyr::select(-Mean, -SD)|>
      pivot_wider(names_from = names, values_from = point)|> 
      dplyr::select(-num)
  }
  if (z_score){
    df_tot = df_tot  |>  
      mutate(across(all_of(labs), scale))
  }
  
  return(df_tot)
  
}

# generate hvs ----
# load all data
d = read_csv('data/CESImixResults.csv')
d
# number or iterations
reps = 3

# generate points and z-score across iterations
set.seed(14)
df = d |> 
  # duplicate for number of reps
  slice(rep(1:n(), each=reps))|> 
  mutate(i = rep(1:reps, times=nrow(d))) |> 
  group_by(i) |> 
  nest() |> 
  # apply function to generate random points
  mutate(points = map(data, \(data) HVvalues(df = data, ID_rows = c('species', 'season'),
                                             names = c('source'), 
                                             mean = 'mean', sd = 'sd', n = 30,
                                             end_points = T, low = 'lowend', high = 'highend',
                                             z_score = T))) |> 
  select(i, points) |> 
  unnest(points)


df

# generate hypervolumes
df = df |> 
  group_by(species, season, i) |> 
  nest() |> 
  mutate(hv = map(data, \(data) hypervolume_gaussian(data, name = paste(species, season, i,sep = '_'),
                                                     samples.per.point = 500,
                                                     kde.bandwidth = estimate_bandwidth(data), 
                                                     sd.count = 3, 
                                                     quantile.requested = 0.95, 
                                                     quantile.requested.type = "probability", 
                                                     chunk.size = 1000, 
                                                     verbose = F)),
         hv_size = map_dbl(hv, \(hv) get_volume(hv)),
         centroid = map(hv, \(hv) get_centroid(hv)))

# df1 = df
#write_rds(df, 'data/CESIhvAll.rds', compress = 'gz')

# overlap within species by season 
ov_sn = df1 |> 
  select(species, season, hv, hv_size, i) |> 
  pivot_wider(names_from = season, values_from = c(hv,hv_size)) |> 
  mutate(size_rat = hv_size_Dry/hv_size_Wet,
         set = map2(hv_Wet,hv_Dry, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv_Wet,hv_Dry, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F))) |> 
  unnest_wider(ov) |> 
  select(species, i, hv_size_Wet, hv_size_Dry, 
         size_rat, jaccard, sorensen,
         uniq_Wet = frac_unique_1, uniq_Dry = frac_unique_2, 
         dist_cent)
ov_sn
#write_csv(ov_sn, 'data/hvOv_season.csv')

#ov_sn = read_csv('data/hvOv_season.csv')

# overlap
cols = c("Pinfish" = 'yellow2',
         "Mojarra" = 'slategray4',
         "Silver perch" = 'snow3',
         "Bay anchovy" = 'deepskyblue1',
         "Pigfish" = 'orange', 
         "Pink shrimp" = 'pink',
         "Rainwater killifish" = 'firebrick',
         'All' = 'black')

df = ov_sn|> 
  group_by(species) |> 
  summarize(mean = mean(sorensen),
            median = median(sorensen),
            low = quantile(sorensen, 0.025),
            up = quantile(sorensen, 0.975))

ggplot(df, aes(x = species, y = mean, color = species))+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2, 
                  position=position_dodge(width = 0.5))+
  labs(x = NULL, y = 'Niche overlap') +
  scale_x_discrete(labels = c("Bay \nanchovy",
                              "Mojarra",
                              "Pigfish",
                              "Pinfish",
                              "Pink \nshrimp",
                              "Rainwater \nkillifish",
                              "Silver \nperch" ))+
  scale_color_manual(values = cols)+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# ggsave("ovPres.png", units="in", width=7, height=6, dpi=600)

# percent unique 
df = ov_sn|>  
  pivot_longer(uniq_Wet:uniq_Dry,names_to = 'season',
               values_to = 'vol') |> 
  group_by(species, season) |> 
  summarize(mean = mean(vol),
            median = median(vol),
            low = quantile(vol, 0.025),
            up = quantile(vol, 0.975)) |> 
  mutate(season = factor(season, levels = c('uniq_Wet', 
                                            'uniq_Dry')))

ggplot(df, aes(x = species, y = mean, color = season))+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2, 
                  position=position_dodge(width = 0.5))+
  labs(x = 'Species', y = 'Niche volume unique',
       color = 'Season') +
  scale_fill_manual(values = c('uniq_Wet' = 'skyblue3', 
                               'uniq_Dry' = 'indianred3'),
                    labels = c('uniq_Wet' = 'Wet', 
                               'uniq_Dry' = 'Dry')) +
  scale_color_manual(values = c('uniq_Wet' = 'skyblue3', 
                                'uniq_Dry' = 'indianred3'),
                     labels = c('uniq_Wet' = 'Wet',
                                'uniq_Dry' = 'Dry')) +
  scale_x_discrete(labels = c("Bay \nanchovy",
                              "Mojarra",
                              "Pigfish",
                              "Pinfish",
                              "Pink \nshrimp",
                              "Rainwater \nkillifish",
                              "Silver \nperch" ))+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# ggsave("uniquePres.png", units="in", width=10, height=7, dpi=600)


# plot size
df = read_csv('data/hvOv_season.csv') |> 
  pivot_longer(hv_size_Wet:hv_size_Dry,names_to = 'season',
               values_to = 'vol') |> 
  group_by(species, season) |> 
  summarize(mean = mean(vol),
            median = median(vol),
            low = quantile(vol, 0.025),
            up = quantile(vol, 0.975)) |> 
  mutate(season = factor(season, levels = c('hv_size_Wet', 
                                            'hv_size_Dry')))

ggplot(df, aes(x = species, y = mean, color = season))+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2, 
                  position=position_dodge(width = 0.5))+
  labs(x = 'Species', y = 'Trophic niche width',
       color = 'Season') +
  scale_fill_manual(values = c('hv_size_Wet' = 'skyblue3', 
                               'hv_size_Dry' = 'indianred3'),
                    labels = c('hv_size_Wet' = 'Wet', 
                               'hv_size_Dry' = 'Dry')) +
  scale_color_manual(values = c('hv_size_Wet' = 'skyblue3', 
                                'hv_size_Dry' = 'indianred3'),
                     labels = c('hv_size_Wet' = 'Wet',
                                'hv_size_Dry' = 'Dry')) +
  scale_x_discrete(labels = c("Bay \nanchovy",
                              "Mojarra",
                              "Pigfish",
                              "Pinfish",
                              "Pink \nshrimp",
                              "Rainwater \nkillifish",
                              "Silver \nperch" ))+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))
# 
# ggsave("sizePres.png", units="in", width=10, height=7, dpi=600)


# centroid distance
df = ov_sn|> 
  group_by(species) |> 
  summarize(mean = mean(dist_cent),
            median = median(dist_cent),
            low = quantile(dist_cent, 0.025),
            up = quantile(dist_cent, 0.975))

ggplot(df, aes(x = species, y = mean, color = species))+
  geom_hline(aes(yintercept = 1), linewidth = 1, linetype = 'dashed')+
  geom_pointrange(aes(ymin = low, ymax = up),
                  size = 1.5, linewidth = 1.5, fatten = 2, 
                  position=position_dodge(width = 0.5))+
  labs(x = NULL, y = 'Centroid distance') +
  scale_x_discrete(labels = c("Bay \nanchovy",
                              "Mojarra",
                              "Pigfish",
                              "Pinfish",
                              "Pink \nshrimp",
                              "Rainwater \nkillifish",
                              "Silver \nperch" ))+
  scale_color_manual(values = cols)+
  scale_y_continuous(limits = c(0,4.1))+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

# ggsave("cdPres.png", units="in", width=7, height=6, dpi=600)


