Testing differences in gene abundance among regions and pathotypes
================

# Introduction

The purpose of this analysis is to assess potential differences in the
distribution of a specific gene among microbes from a specific microbe.
Because this document is public, I’m not going to name the gene or the
microbe.

## ANOVA

A simple way to do this would be an ANOVA. Let’s try and check whether
the distributions of residuals and stability of variances look
reasonable.

``` r
library(tidyverse)
library(MASS)
theme_set(theme_classic()) 
library(furrr)
library(tictoc)
library(rlang) # this is maybe needed by calculate_mean_diff
plan(multisession)

# Read in raw zor-orz data
region <- readxl::read_xlsx("data/FINAL data for Steen.xlsx", 
                            sheet = "zor-orz gene number",
                            range = "A11:G18") %>%
  rename(gene.count = Continent) %>%
  pivot_longer(-1, names_to = "category", values_to = "count") %>%
  group_by(category) %>%
  mutate(freq = count / sum(count, na.rm = TRUE))
  
path <- readxl::read_xlsx("data/FINAL data for Steen.xlsx", 
                            sheet = "zor-orz gene number",
                            range = "A1:E8") %>%
  rename(gene.count = Pathotype) %>%
  pivot_longer(-1, names_to = "category", values_to = "count") %>%
  group_by(category) %>%
  mutate(freq = count / sum(count, na.rm = TRUE))


# path <- read_csv("data/pathotype.csv") %>%
#   #filter(`gene number`!="TOTAL") %>%
#   pivot_longer(cols = 2:5, names_to = "category", values_to="count") %>%
#   mutate(gene.count = as.numeric(`gene number`)) %>%
#   dplyr::select(-`gene number`) %>%
#   group_by(category) %>%
#   mutate(freq = count / sum(count, na.rm = TRUE))

# I need to recreate the raw data to do an anova
# I wrote some *very* ugly code to do this, so I hid it in a separate file
source("R/helper_funs.R")


raw_region_data <- recreate_raw(region) %>%
  arrange(category) # this appears to have worked
raw_path_data <- recreate_raw(path)

raw_plot <- function(df) {
  p <- ggplot(df, aes(x=category, y=gene.count)) + 
  geom_boxplot() + 
  #geom_point(position=position_jitter(height = 0.3), alpha = 0.5) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle=-45, hjust=0))
  p
}

p_region <- raw_plot(raw_region_data)
p_path <- raw_plot(raw_path_data)

raw_line_plot <- function(df) {
  p <- ggplot(df, aes(x=gene.count, y=freq, color=category)) + 
    geom_point() + 
    geom_line()  + 
    scale_color_brewer(type="qual")
  p
}

raw_bar_plot <- function(df) {
  p <- ggplot(df, aes(x=gene.count, y=freq, fill=category)) + 
    geom_bar(stat="identity", position="dodge") + 
    scale_fill_brewer(type="qual")
  p
}

# Make up plots for later display
p_reg_line <- raw_line_plot(raw_region_data)
p_path_line <- raw_line_plot(raw_path_data)
p_reg_bar <- raw_bar_plot(raw_region_data)
p_path_bar <- raw_bar_plot(raw_path_data)
```

I’m not really sure what the best way to display these is, so I’m giving
three options:

![](selene_stats_files/figure-commonmark/unnamed-chunk-2-1.png)

These data sets have some subtle differences in distributions among
pathotypes/regions, which are easiest to see in the middle row of
colored line plots. Is a linear model (ANOVA) good for these data?
Specifically: ANOVA is fairly robust to unbalanced designs and to
heteroskedasticity, but not to hetereoskedastic data in an unbalanced
design. So let’s check the heteroskedasticity.

### Linear model for regions

``` r
region_model <- lm(gene.count ~ category, data=raw_region_data)
summary(aov(region_model))
```

                  Df Sum Sq Mean Sq F value   Pr(>F)    
    category       5     90  17.983   9.335 8.61e-09 ***
    Residuals   1780   3429   1.926                     
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

This model finds signficant differences among regions. But before we
take this too seriously, let’s check whether the residuals are normally
distributed. A good way to do that is via a QQ plot. The

``` r
plot(region_model, which=2) 
```

![](selene_stats_files/figure-commonmark/unnamed-chunk-4-1.png)

I’d say we these residuals are sufficiently non-normally distributed
that a linear model is not a good choice for the regional data.

![](https://64.media.tumblr.com/f8b7a14c2fa304a712b5f92ea14d62f9/tumblr_n41bxrhleZ1rvirvyo1_400.gif)

### Linear model for pathotypes

``` r
path_model <- lm(gene.count ~ category, data=raw_path_data)
summary(aov(path_model))
```

                 Df Sum Sq Mean Sq F value Pr(>F)    
    category      3  201.9   67.30   46.03 <2e-16 ***
    Residuals   601  878.6    1.46                   
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Again, significant differences among pathotypes.

``` r
plot(path_model, which=2)
```

![](selene_stats_files/figure-commonmark/unnamed-chunk-6-1.png)

Same situation here. The QQ-plot is sufficiently
not-like-a-straight-line that I don’t really want to interpret the p
values that come from it.

## Poisson distribution?

I propose that we can think of gene distribution as a poisson process,
where different values of $\lambda$ indicate different probabilities of
the gene being “handed out”. If this is the case, we can assess whether
there are difference in lambda among regions or pathotypes - but first
we need to assess whether the data are, in fact, poisson-distributed.
We’ll simply load the data, fit it to a poisson distribution, and see
whether the fit looks good. I think in this case a qualitative
assessment is at least as good as some kind fo statistical test of
goodness-of-fit.

``` r
d <- rbind(region %>% mutate(type="region"),
           path %>% mutate(type="path")) %>%
  mutate(category=factor(category, levels = c("Africa","Asia", "Europe", "North America", "Oceania", "South America", "pa", "healthy", "intestinal disease", "urinary disease"), ordere=TRUE))

ggplot(d, aes(x=gene.count, y=count)) + 
  geom_point() + 
  facet_wrap(~category, scale="free_y")
```

![](selene_stats_files/figure-commonmark/unnamed-chunk-7-1.png)

These data do not look poisson-distributed. I’m pretty sure that part of
the issue is there is correlation between the two genes in terms of
whether they are likely to appear in hte genome - that is, if one of the
genes is present, the other is likely to be as well. Note there are
almost never exactly 3 genes present. So, I don’t really want to model
this with Poisson distributions.

# Monte Carlo simulations

I think it is more robust to do a Monte Carlo simulation of variation in
the ANOVA *f* ratio.

``` r
# There are, like, a lot of faster ways to do this

n <- 1000
nrow.reg <- nrow(region)
nrow.path <- nrow(path)
# reg.f.vec <- vector("double", n)
# path.f.vec <- vector("double", n)

#set.seed(512)
# region.loop.time <- system.time({
#   for(i in 1:n) {
#   reg.f.vec[i] <- shuf_calc_f(region, nrow.reg)
# }
# })

tic()
reg.f.vec <- future_map_dbl(seq_along(1:n), shuf_calc_f, df=region, nrow=nrow.reg, 
                             .options = furrr_options(seed = TRUE)) 
path.f.vec <- future_map_dbl(seq_along(1:n), shuf_calc_f, df=path, nrow=nrow.path, 
                             .options = furrr_options(seed = TRUE)) 
toc()
```

    2.882 sec elapsed

``` r
# OK, I have to ask: is this faster than map_dbl? THe time to beat is 12.37 (I know, I shoudl repeat it...)
tic()
reg.f.vec <- map_dbl(seq_along(1:n), shuf_calc_f, df=region, nrow=nrow.reg, 
                             .options = furrr_options(seed = TRUE)) 
path.f.vec <- map_dbl(seq_along(1:n), shuf_calc_f, df=path, nrow=nrow.path, 
                             .options = furrr_options(seed = TRUE)) 
toc() # Yep, seems like 4 times faster! Which is good, because it is trying to use 6 cores, and 6 = 4 (kind of )
```

    5.608 sec elapsed

``` r
# tic.clearlog()
# tic()
# path.loop.time <- system.time({
#   for(i in 1:n) {
#     path.f.vec[i] <- shuf_calc_f(path, nrow.path)
#   }
# })
# loop.time <- toc() # serial loop is maybe ~4 times slower than parallel method


# This takes ~22 sec per loop on my system

f_vals <- data.frame(reg.sim.f = reg.f.vec,
                     path.sim.f = path.f.vec)
```

``` r
# Pull out actual f values
reg.f.real <- summary(aov(region_model))[[1]][1,4]
path.f.real <- summary(aov(path_model))[[1]][1,4]
```

How do the real f values compare to the simulated, null-hypothesis
values?

``` r
p_reg_hist <- ggplot(f_vals, aes(x=reg.sim.f)) + 
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = reg.f.real, color="red")  + 
  ggtitle("data by region")
print(p_reg_hist)
```

![](selene_stats_files/figure-commonmark/unnamed-chunk-11-1.png)

``` r
p_path_hist <- ggplot(f_vals, aes(x=path.sim.f)) + 
  geom_histogram(bins = 100) + 
  geom_vline(xintercept = path.f.real, color="red") + 
  ggtitle("data by pathology")
print(p_path_hist)
```

![](selene_stats_files/figure-commonmark/unnamed-chunk-12-1.png)

So: I have simulated 1,000 and found that, for each case, the actual
measured *f* values are much, much larger than they would be likely to
be if the null hypothesis were true - so much larger that we can’t
calculate a p value, because none of our 1,000 simulations captured a
*f* value that big. We can say conservatively say that, in each case, p
\< 0.0001.

In summary:

| f-value      | Simulated maximum | Observed |
|--------------|-------------------|----------|
| by pathology | 8.1               | 46       |
| by region    | 5.2               | 9.3      |

# Monte Carlo Tukey Test

A Tukey test works by comparing means of all possible combinations of
populations (in this case, regions or pathotypes) and then comparing to
a studentized range distribution. I’m going to do exactly this, except
that the studentized range distribution is replaced with the observed
distribution of mean differences in shuffled data.

``` r
# mean_diff <- function(x) { # analagous to 
#   abs(mean(x[[1]], na.rm = TRUE) - mean(x[[2]], na.rm = TRUE))
# }

# path_means <- raw_path_data %>% # I feel like this is duplicative but I guess not?
#   # I know I could calulate it from path, but it feels like it is too easy to screw that up even though it requires roughly 6th grade math
#   group_by(category) %>%
#   summarise(mean.gene.count = mean(gene.count, na.rm=TRUE))




actual_path_mean_diffs <- shuf_and_calc_means(raw_path_data, shuf = FALSE) # Seems to work
n <- 10
null_path_mean_diffs <- future_map(seq_along(1:n), 
                                   shuf_and_calc_means, df=raw_path_data, shuf=TRUE,
                                   .options = furrr_options(seed = TRUE)) 

# This takes the list from null_path_mean_list and turns it into a tidy data set
null_path_distribs <- sim_list_to_df(null_path_mean_diffs, summarise = TRUE, alpha = 0.95)
merge(null_path_distribs, actual_path_mean_diffs, by = "dif.id") %>%
  mutate(sig.diff = case_when(mean.diff >= cutoff.diff ~ TRUE, TRUE ~ FALSE))
```

                                  dif.id cutoff.diff  mean.diff sig.diff
    1                 bacteremia-healthy   0.4906021 0.26840634    FALSE
    2      bacteremia-intestinal disease   0.3444305 1.34124149     TRUE
    3         bacteremia-urinary disease   0.4871658 0.09124721    FALSE
    4         healthy-intestinal disease   0.3388270 1.07283515     TRUE
    5            healthy-urinary disease   0.3867725 0.17715913    FALSE
    6 intestinal disease-urinary disease   0.2387823 1.24999428     TRUE

Now I write a function to shuffle the data and perform the same
calculation, and then return min and max values. (Actually I guess i
just need the max value?)

More accurately the approach here is to:

- For each pair of comparisons, compare the actual mean difference to
  the set of mean differences for shuffled data. There’s no need to
  correct for standard error (I don’t think).

``` r
# shuf_and_calc_means <- function(df, shuf = FALSE) {
#   # Shuffle data 
#   if(shuf) {
#     df$category <- sample(df$category, size = nrow(df), replace = TRUE)
#   }
#   
#   
#   # Calculate means
#   means <- df %>%
#     group_by(category) %>%
#     summarise(mean.gene.count = mean(gene.count, na.rm = TRUE))
#   
#   # Calculate mean differences for each set of groups
#   diffs <- calc_mean_diff(means)
# }
```
