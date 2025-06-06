# title: "Useful functions for hypervolumes"

# 
# This script contains information about the functions used in the hypervolume scripts.  
# 

## Merge/Join
# If two data frames contain different columns of data, then they can be merged together with the family of join functions.
# 
# +`left_join()` = uses left df as template and joins all matching columns from right df 
# +`right_join()` = uses right df as template and joins all matching columns from left df
# +`inner_join()` = only matches columns contained in both dfs
# +`full_join()` = combines all rows in both dfs

library(tidyverse)

left = tibble(name = c('a', 'b', 'c'),
              n = c(1, 6, 7), 
              bio = c(100, 43, 57))

right = tibble(name = c('a', 'b', 'd', 'e'),
               cals = c(500, 450, 570, 600))

left_join(left, right, by = 'name')

right_join(left, right, by = 'name')

inner_join(left, right, by = 'name')

full_join(left, right, by = 'name')

# multiple matches
fish = tibble(species = rep(c('Salmon', 'Cod'),times = 3),
              year = rep(c(1999,2005,2020), each = 2),
              catch = c(50, 60, 40, 50, 60, 100))

col = tibble(species = c('Salmon', 'Cod'),
             coast = c('West', 'East'))

left_join(fish, col, by = 'species')


# ## scaling data
# Because hypervolumes can be generated with any continuous data as an axes, many of the times the units are not comparable. Blonder et al. [2014](https://doi-org.ezproxy.fiu.edu/10.1111/geb.12146) & [2018](https://doi-org.ezproxy.fiu.edu/10.1111/2041-210X.12865) to convert all of the axes into the same units. This can be done by taking the z-score of the values to convert units into standard deviations. Z-scoring data can be done with the formula:
#   $$ z = \frac{x_{i}-\overline{x}}{sd} $$ Where $x_{i}$ is a value, $\overline{x}$ is the mean, and $sd$ is the standard deviation. By z-scoring each axis, 0 is the mean of that axis, a value of 1 means that the value is 1 standard deviation above the global mean of that axis, and a value of -1 is 1 standard deviation below the global mean of the axis. In R this can be done manually or with the `scale()` function. 


fish = tibble(species = rep(c('Salmon', 'Cod'),times = 3),
              year = rep(c(1999,2005,2020), each = 2),
              catch = c(50, 60, 40, 50, 60, 100))

#
fish = fish |> 
  mutate(zcatch1 = (catch - mean(catch))/sd(catch), # manual
         zcatch2 = scale(catch)) # with scale

fish |> group_by(species) |> 
  mutate(zcatch1 = (catch - mean(catch))/sd(catch), # manual
         zcatch2 = scale(catch))
fish 

# center = mean, scale = sd
fish$zcatch2


## nesting data
# One benefit of `tibbles` is that they can contain list columns. This means that we can make columns of `tibbles` that are nested within a dataset. Nesting creates a list-column of data frames; unnesting flattens it back out into regular columns. Nesting is a implicitly summarising operation: you get one row for each group defined by the non-nested columns. This is useful in conjunction with other summaries that work with whole datasets, most notably models. This can be done with the `nest()` and then flattened with `unnest()`

fish = tibble(species = rep(c('Salmon', 'Cod'),times = 3),
              year = rep(c(1999,2005,2020), each = 2),
              catch = c(50, 60, 40, 50, 60, 100))

# using group_by
fish_nest = fish |> 
  group_by(species) |> 
  nest()

fish_nest
fish_nest$data
fish_nest$data[[1]]

# using .by in nest
# column name becomes data unless you change .key
fish_nest2 = fish |> 
  nest(.by = year, .key = 'df')

fish_nest2
fish_nest2$df

# ## map
# ### `purr`
# The newest and new standard package with `tidyverse` is `purr` with its set of `map()` functions. Some similarity to `plyr` (and base) and `dplyr` functions but with more consistent names and arguments. Notice that map function can have some specification for the type of output.
# + `map()` makes a list.
# + `map_lgl()` makes a logical vector.
# + `map_int()` makes an integer vector.
# + `map_dbl()` makes a double vector.
# + `map_chr()` makes a character vector.

df = iris  |> 
  select(-Species)
#summary statistics
map_dbl(df, mean)

# using map with mutate and nest
d = tibble(species = rep(c('Salmon', 'Cod'),times = 3),
           year = rep(c(1999,2005,2020), each = 2),
           catch = c(50, 60, 40, 50, 60, 100)) |> 
  nest(.by = species) |> 
  mutate(correlation = map(data, \(data) cor.test(data$year, data$catch)))

d
d$correlation

# extract information from list
d = d |> 
  mutate(r = map_dbl(correlation, \(x) x$estimate),
         p = map_dbl(correlation, \(x) x$p.value))

d
# Sometimes, there are multiple arguments required for the function you are mapping. You can provide two lists (or columns) using `map2()`. If you have more than two lists needed, you can use `pmap()`. These functions work the same as `map()` based on the desired output. 

library(performance)
# using map2 and pmap
df = iris  |> 
  nest(.by = Species) |> 
  mutate(sw = map(data, \(data) lm(Sepal.Length ~ Sepal.Width, data = data)),
         pl = map(data, \(data) lm(Sepal.Length ~ Petal.Length, data = data)),
         mod_c = map2(sw,pl, \(sw,pl) compare_performance(sw,pl)),
         top = map_chr(mod_c, \(x) x$Name[which.min(x$AICc)]),
         sum_top = pmap(list(sw,pl,top), \(x,y,z) 
                        if (z == 'sw'){
                          summary(x)
                        }else{
                          summary(y)}),
         p = map_dbl(sum_top, \(s) s$coefficients[2, 4]))


df

