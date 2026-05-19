# Libraries

### *Data Manipulation*

```{r}
library(tidyverse)
library(stringr)
library(purrr)
```

### *Analysis*

```{r}
library(psych)
library(ReX)

library(lme4)
library(lmerTest)
library(performance)

## See Hehrman & Xie (2021)|https://osf.io/anwx2/
source("bootstrap_icc.R")
```

### *Plotting*

```{r}
library(ggplot2)
library(sjPlot)
library(tidytext)

library(RColorBrewer)
library(viridis)
library(hrbrthemes)
library(pals)
```

# Functions

### *plotting*

```{r}
## Removes leading zeros
scaleFUN <- function(x) sprintf("%.2f", x)

## Standard Error
se <- function(x) {
  
  n <- length(x)
  stdev <- sd(x)
  
  return(stdev/sqrt(n))
  
}

## 95% Confidence Intervals
sample_lCI <- function(x) {
  
  xbar <- mean(x)
  n <- length(x)
  se <- se(x)
  
  tcrit <- qt(0.975,df=n-1)
  margin <- tcrit*se
  
  return(xbar - margin)
}

sample_uCI <- function(x) {
  
  xbar <- mean(x)
  n <- length(x)
  se <- se(x)
  
  tcrit <- qt(0.975,df=n-1)
  margin <- tcrit*se
  
  return(xbar + margin)
}

## Percentile-based bootstrapped 95% CIs
boot_uCI <- function(x){
  tCI <- numeric(2)
  x_bar <- mean(x)
  
  # lCI <- unname(2 * x_bar - quantile(x, c(0.975, 0.025)))[1]
  uCI <- unname(2 * x_bar - quantile(x, c(0.975, 0.025)))[2]
  
  return(uCI)
}

boot_lCI <- function(x){
  tCI <- numeric(2)
  x_bar <- mean(x)
  
  lCI <- unname(2 * x_bar - quantile(x, c(0.975, 0.025)))[1]
  # uCI <- unname(2 * x_bar - quantile(x, c(0.975, 0.025)))[2]
  
  return(lCI)
}
```

# Dataset

```{r}
master_set <- read.csv('master_set.csv')
```

## Parallel Analysis

```{r}
test_pa <- fa.parallel(select(original, contains("_response")), fa = "pc", sim = F)


scree_data <- as.data.frame(test_pa$pc.values) %>%
  mutate(pc_number = row.names(.)) %>%
#  cbind(as.data.frame(test_pa$pc.sim)) %>%
  cbind(as.data.frame(test_pa$pc.simr)) %>%
  `colnames<-`(c("real_vals", "compnum", "resamp_vals")) %>%
  .[,c(2, 1, 3)] %>%
  pivot_longer(real_vals:resamp_vals,
                names_to = "pc_type",
                values_to = "eigenvalue")

scree_data$exp_var <- (scree_data$real_vals/16)*100

paste("explained variance =", sum(scree_data[1:4,'exp_var']))
```

## Scree Plot

```{r}
scree_plot <- ggplot(scree_data, aes(x = as.numeric(compnum), y = eigenvalue, group = pc_type, color = pc_type, linetype = pc_type)) +
  geom_point(show.legend = T, size = 2) +
  geom_line(size = 1) +
  labs(x = "Component Number",
       y = "Eigenvalue") +
  scale_color_manual(values = c("dodgerblue", "red2"),
                     name = "",
                     labels = c("Real Data", "Resampled Data")) +
  scale_linetype_discrete(name = "",
                     labels = c("Real Data", "Resampled Data")) +
  scale_x_continuous(breaks = seq(1,16)) +
  theme_sjplot() +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10, color = "black"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black"),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold"))

scree_plot
```

## Task-level mean component scores

```{r}
## Colour-coding task categories (uses RColorBrewer)
category_palettes <- list(
  ExecutiveWM = "Oranges",
  Visuomotor = "Greens",
  InternallyOriented = "Blues",
  Naturalistic = "Purples"
)

category_base_colors <- map_chr(category_palettes, ~ brewer.pal(5, .x)[3])
```

### *Wrangling*

```{r}
task_avs <- master_set %>%
  select(., c('commspsych', 'PC1', 'PC2', 'PC3', 'PC4')) %>%
  group_by(commspsych, .drop = F) %>%
  summarise(across(PC1:PC4, mean, .names = "mean_{.col}"),
            across(PC1:PC4, sd, .names = "sd_{.col}"),
            across(PC1:PC4, sample_lCI, .names = "lCI_{.col}"),
            across(PC1:PC4, sample_uCI, .names = "uCI_{.col}"),
            across(PC1:PC4, se, .names = "se_{.col}")) %>%
  mutate(task_cat = case_when(
    str_detect(commspsych, "^2B|Math$|1B") ~ "ExecutiveWM",
    str_detect(commspsych, "0B|FingerTap|GoNoGo") ~ "Visuomotor",
    str_detect(commspsych, "Read|Documentary|SciFi") ~ "Naturalistic",
    str_detect(commspsych, "Memory|Friend|You") ~ "InternallyOriented"
  )) %>%
  transform(commspsych=factor(commspsych, levels=c(
    "HardMath","EasyMath",
    "2B-Scene","2B-Face","1B",
    "0B","GoNoGo","FingerTap",
    "You","Friend","Memory",
    "Read","Documentary","SciFi"))) %>%
  arrange(commspsych) %>%
  mutate(taskord = seq(levels(commspsych)))


task_avs_colored <- task_avs %>%
  group_by(task_cat) %>%
  mutate(n_tasks = n()) %>%
  mutate(color = brewer.pal(max(n_tasks), category_palettes[[unique(task_cat)]]),
         category_fill = category_base_colors[[unique(task_cat)]]) %>%
  ungroup() %>%
  pivot_longer(mean_PC1:se_PC4,
               names_to = c(".value", "comp"),
               names_sep="_" ) %>%
  mutate(., label = case_when(
    comp == 'PC1' ~ 'Episodic Knowledge',
    comp == 'PC2' ~ 'Intrusive Distraction',
    comp == 'PC3' ~ 'Deliberate Task-Focus',
    comp == 'PC4' ~ 'Sensory Engagement'
  )) %>%
  transform(label=factor(label,levels=c(
    "Episodic Knowledge",
    "Intrusive Distraction",
    "Deliberate Task-Focus",
    "Sensory Engagement"))) %>%
  transform(task_cat=factor(task_cat, levels=c(
    "ExecutiveWM",
    "Visuomotor",
    "Naturalistic",
    "InternallyOriented"
  ))) %>%
  mutate(catord = case_when(
    task_cat == "ExecutiveWM" ~ 1,
    task_cat == "Visuomotor" ~ 2,
    task_cat == "Naturalistic" ~ 4,
    task_cat == "InternallyOriented" ~ 3
  ))
```

### *Plotting*

```{r}
ggplot(task_avs_colored,
       aes(x = commspsych,
           y = scale(mean, scale = F),
           fill = color,
           color = task_cat)) +
  geom_col() +
  geom_errorbar(aes(ymin = lCI,
                    ymax = uCI),
                width = 0.5, size = 0.95, color = "black") +
  facet_grid(~ label, scales = "free_x") +
  scale_fill_identity() + 
  ylab("Component Score") +
  xlab("Task") +
  ylim(-1.75, 1.75) +
  coord_flip() +
  theme_sjplot() +
  theme(axis.text.x = element_text(face = "bold", color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold", color = "black",
                                   size = 12),
        axis.title.x = element_text(face = "bold", size = 12, color = "black"),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(face = "bold", color = "black", size = 12),
        strip.text = element_blank()) +
  scale_color_manual("",
                     limits=unique(task_avs_colored$task_cat),
                     labels=c("Executive/WM","Visuomotor",
                              "Internally-Oriented","Naturalistic"),
                     values=unique(task_avs_colored$category_fill),
                     guide = guide_legend(override.aes = list(
                       fill = category_base_colors), ncol = 2))
```

## Task-level mean gradient coordinates

### *Wrangling

```{r}
grad_avs <- master_set %>%
  select(., c('commspsych', contains("gradient"))) %>%
  group_by(commspsych, .drop = F) %>%
  summarize(across(gradient1:gradient5, mean)) %>%
  mutate(task_cat = case_when(
    str_detect(commspsych, "^2B|Math$|1B") ~ "ExecutiveWM",
    str_detect(commspsych, "0B|FingerTap|GoNoGo") ~ "Visuomotor",
    str_detect(commspsych, "Read|Documentary|SciFi") ~ "Naturalistic",
    str_detect(commspsych, "Memory|Friend|You") ~ "InternallyOriented"
  )) %>%
  transform(commspsych=factor(commspsych, levels=c(
    "HardMath","EasyMath",
    "2B-Scene","2B-Face","1B",
    "0B","GoNoGo","FingerTap",
    "You","Friend","Memory",
    "Read","Documentary","SciFi"))) %>%
  arrange(commspsych) %>%
  mutate(taskord = seq(levels(commspsych)))

grad_avs_colored <- grad_avs %>%
  group_by(task_cat) %>%
  mutate(n_tasks = n()) %>%
  mutate(color = brewer.pal(max(n_tasks), category_palettes[[unique(task_cat)]]),
         category_fill = category_base_colors[[unique(task_cat)]]) %>%
  ungroup() %>%
  pivot_longer(gradient1:gradient5, names_to = 'gradient', values_to = 'coord') %>%
  mutate(label = case_when(
    gradient == 'gradient1' ~ 'Sens. – Assoc.',
    gradient == 'gradient2' ~ 'Mot. – Vis.',
    gradient == 'gradient3' ~ 'DMN – Cont.',
    gradient == 'gradient4' ~ 'DAN – Sens.Att.',
    gradient == 'gradient5' ~ 'Sens.Mot. – Verb.Sem.'
  )) %>%
  transform(label=factor(label,levels=c(
    'Sens. – Assoc.',
    'Mot. – Vis.',
    'DMN – Cont.',
    'DAN – Sens.Att.',
    'Sens.Mot. – Verb.Sem.'
    )))
```

### *Plotting*

```{r}
ggplot(grad_avs_colored,
       aes(x = Task_name,
           y = scale(coord, scale = F),
           fill = color,
           color = task_cat)) +
  geom_col() +
  facet_wrap(~ label, scales = "free_x", ncol = 5, nrow = 1) +
  scale_fill_identity() + 
  ylab("Coordinate") +
  xlab("Task") +
  # scale_x_reordered() +
  # ylim(-0.5, 0.5) +
  scale_y_continuous(limits = c(-0.4, 0.4),
                     breaks = c(-0.25, 0, 0.25)) +
  coord_flip() +
  theme_sjplot() +
  theme(axis.text.x = element_text(face = "bold", color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold", color = "black",
                                   size = 12),
        axis.title.x = element_text(face = "bold", size = 12, color = "black"),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.title = element_blank(),
        strip.background = element_rect(fill=NA,colour=NA),
        strip.text = element_text(face="bold", color = "black", size = 12)) +
  scale_color_manual("Task Category",
                     limits=unique(task_avs_colored$task_cat),
                     labels=c("Executive/WM","Perceptual",
                              "Internally-Oriented","Naturalistic"),
                     values=unique(task_avs_colored$category_fill),
                     guide = guide_legend(override.aes = list(fill = category_base_colors)))
```


## Bootstrapped ICC values

### Packaging Hehrman & Xie (2021)'s bootICC() function

```{r}

model = NA
iterations = 1000
boot.type = "parametric"
boot.ci.type = "basic"
reps = F

run_bootICC <- function(dv){
  
  na_inds <- which(is.na(dv))
  model_df <- formdf[!(row.names(formdf) %in% na_inds), ]
    
  dv <- dv[!is.na(dv)]
  
  
  model <- lmer(as.formula(paste(formula, collapse = " ")), model_df)
  bootICC <- bootstrapICC(model, iterations, boot.type, boot.ci.type)
  bootframe <- as.data.frame(bootICC$CIs)
  
  if (reps) {
      repframe <- as.data.frame(bootICC$boot_object$replicates)
      frames = list("iccs" = bootframe, "reps" = repframe)
      return(frames)
  }
  
  return(bootframe)
}

```

### Mean task-level component scores per subject

```{r}
subject_avgs <- master_set %>% 
  group_by(subID, commspsych, .drop = F) %>%
  summarise(PC1 = mean(PC1),
            PC2 = mean(PC2),
            PC3 = mean(PC3),
            PC4 = mean(PC4))
```

### Trait-like (Subject-level) stability

```{r}

model = 'dv ~ 1 + (1|subID)'
model_df = subject_avgs
reps = F

traitICCs <- do.call("rbind", apply(subject_avgs[3:6], 2, run_bootICC))

## To save replicates
# names(traitICC_reps) <- c('PC1', 'PC2', 'PC3', 'PC4')
# traitreps <- map_dfr(names(traitICC_reps), function(pc) {
#   reps <- traitICC_reps[[pc]]
# 
#   cbind(data.frame(component = pc), reps)
#   
# })

traitICCs_wide <- traitICCs %>%
  mutate(indices = row.names(.)) %>%
  separate(indices, c("comp", "measure"), sep = "\\.") %>%
  select(., -measure) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  select(., c('comp', 'subID ICC', 'lower', 'upper')) %>%
  `colnames<-`(c("thought_pattern", "subject_ICC", "LCI", "UCI")) %>%
  drop_na() %>%
  mutate(thought_pattern = c('Episodic Knowledge', 'Intrusive Distraction', 'Deliberate Task-Focus', 'Sensory Engagement'))


traitICCs_wide
```

### State-like (task-level) Stability

```{r}

model = 'dv ~ 1 + (1|commspsych)'
model_df = subject_avgs
reps = F
 
taskwiseICCs <- do.call("rbind", apply(subject_avgs[3:6], 2, run_bootICC))

taskwiseICCs_wide <- taskwiseICCs %>%
  mutate(indices = row.names(.)) %>%
  separate(indices, c("comp", "measure"), sep = "\\.") %>%
  select(., -measure) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  select(., c('comp', 'commspsych ICC', 'lower', 'upper')) %>%
  `colnames<-`(c("thought_pattern", "task_ICC", "LCI", "UCI")) %>%
  drop_na() %>%
  mutate(thought_pattern = c('Episodic Knowledge', 'Intrusive Distraction', 'Deliberate Task-Focus', 'Sensory Engagement'))


taskwiseICCs_wide
```

### Cross-factored State x Trait (Subject:Task) Stability

```{r}

idtask <- master_set %>%
  mutate(IDtask = paste(subID, commspsych, sep = ""))

model = 'dv ~ 1 + (1|subID:commspsych)'
model_df = master_set
reps = F

statetraitICCs <- do.call("rbind", apply(master_set[3:6], 2, run_bootICC))

statetraitICCs_wide <- statetraitICCs %>%
  mutate(indices = row.names(.)) %>%
  separate(indices, c("comp", "measure"), sep = "\\.") %>%
  select(., -measure) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  select(., c('comp', 'subID:commspsych ICC', 'lower', 'upper')) %>%
  `colnames<-`(c("thought_pattern", "subject_task_ICC", "LCI", "UCI")) %>%
  drop_na() %>%
  mutate(thought_pattern = c('Episodic Knowledge', 'Intrusive Distraction', 'Deliberate Task-Focus', 'Sensory Engagement'))
  
stateICCs_wide
```

### Trait-like Stability Per Task

```{r}
taskPC <- idtask %>%
  select(., c('subID', 'commspsych', 'IDtask', 'PC1', 'PC2', 'PC3', 'PC4')) %>%
  group_by(IDtask) %>%
  mutate(taskprobe = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = !IDtask,
              names_from = commspsych,
              values_from = c('PC1', 'PC2', 'PC3', 'PC4')) %>%
  `colnames<-`(gsub('-', '_', colnames(.)))


model = 'dv ~ 1 + (1|subID)'
model_df = taskPC
reps = T

pertaskICCs <- do.call("rbind", apply(taskPC[3:ncol(taskPC)], 2, run_bootICC))

write.csv(pertaskICCs, 'bootstrappediccs_task.csv')

pertaskICCs_wide <- pertaskICCs %>%
  mutate(indices = row.names(.)) %>%
  separate(indices, c("taskPC", "measure"), sep = "\\.") %>%
  separate(taskPC, c("PC", "commspsych"), sep = "_", extra = "merge", fill = "left") %>%
  select(., -measure) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  select(., c('PC', "commspsych", 'subID ICC', 'lower', 'upper')) %>%
  `colnames<-`(c("thought_pattern", "task", "subject_ICC", "LCI", "UCI")) %>%
  drop_na() %>%
  mutate(thought_pattern = case_when(
    thought_pattern == 'PC1' ~ 'Episodic Knowledge',
    thought_pattern == 'PC2' ~ 'Intrusive Distraction',
    thought_pattern == 'PC3' ~ 'Deliberate Task-Focus',
    thought_pattern == 'PC4' ~ 'Sensory Engagement'
  ))

write.csv(pertaskICCs_wide, 'bootstrappediccs_task(filtered).csv')
```

#### *Wrangling*

```{r}
pertaskICCs_mtx <- pertaskICCs_wide %>%
  select(., c('thought_pattern', 'task', 'subject_ICC')) %>%
  pivot_wider(names_from = thought_pattern, values_from = subject_ICC) %>%
  `colnames<-`(gsub(' ', '_', colnames(.))) %>%
  mutate(taskmeans = rowMeans(.[2:5]))

taskICCs_colored <- task_avs_colored %>%
  select(., c(commspsych, task_cat, color, category_fill)) %>%
  merge(pertaskICCs_wide, by.x = "commspsych", by.y = "task") %>%
  distinct()
```

#### *Barplots*

```{r}

t_shift <- scales::trans_new("shift",
                             transform = function(x) {x-0.5},
                             inverse = function(x) {x+0.5})

pertaskBars <- ggplot(transform(taskICCs_colored,
                 comp_label=factor(comp_label,
                              levels=c("Episodic Knowledge",
                                           "Intrusive Distraction",
                                           "Deliberate Task-Focus",
                                           "Sensory Engagement"))),
                 aes(x = commspsych,
                           y = estimate, 
                           fill = color,
                           color = task_cat)) +
  geom_col(show.legend = F) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, size = 0.95, color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(y = expression('ICC'[italic('Subject')]~' (Relative to 0.5)'),
       x = "Task") +
  scale_y_continuous(trans = t_shift, limits = c(0.15, 0.85), breaks = c(0.25, 0.50, 0.75)) +
  facet_grid(~ comp_label) +
  scale_fill_identity() +
  # annotate('text',
  #          label=expression(sigma[italic('w')]^{2}~'>'~sigma[italic('b')]^{2}),
  #          x = 7.5, y = 0.1, fontface = 2) +
  # annotate('text',
  #          label=expression(sigma[italic('b')]^{2}~'>'~sigma[italic('w')]^{2}),
           # x = 7.5, y = 0.9, fontface = 2) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(face = 'bold', color = 'black', size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face = 'bold', color = 'black', size = 15),
        axis.text.x = element_text(face = 'bold', color = 'black', size = 9),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_blank(),
        strip.text = element_blank()) +
  scale_color_manual("Task Category",
                     limits=unique(task_avs_colored$task_cat),
                     labels=c("Executive/WM","Perceptual",
                              "Internally-Oriented","Naturalistic"),
                     values=unique(task_avs_colored$category_fill),
                     guide = guide_legend(override.aes = list(fill = category_base_colors)))

pertaskBars

```

#### *Field Maps*

```{r}
taskfield <- pertaskICCs_wide %>%
  subset(., term == "subID Variance" | term == "Residual Variance") %>%
  select(., c('task', 'comp', 'term', 'estimate')) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  mutate(component = case_when(
    comp == 1 ~ 'Episodic Knowledge',
    comp == 2 ~ 'Intrusive Distraction',
    comp == 3 ~ 'Deliberate Task-Focus',
    comp == 4 ~ 'Sensory Engagement'
    )) %>%
  `colnames<-`(c('task', 'PC', 'sigma2_b', 'sigma2_w', 'component'))


taskfield_colored <- taskfield %>%
  merge(labels_master, by.x = "task", by.y = "iccs") %>%
  merge(task_avs_colored, by.x = c("commspsych", "component"), by.y = c("commspsych", "label")) %>%
  distinct() %>%
  transform(component=factor(component, levels=c(
              "Episodic Knowledge",
              "Intrusive Distraction",
              "Deliberate Task-Focus",
              "Sensory Engagement"))) %>%
  transform(task_cat=factor(task_cat, levels=c(
    "ExecutiveWM",
    "Visuomotor",
    "Naturalistic",
    "InternallyOriented"
  )))


pertaskfield <- pertaskICCs %>%
  mutate(indices = row.names(.)) %>%
  separate(indices, c("taskPC", "measure"), sep = "\\.") %>%
  separate(taskPC, c("PC", "task"), sep = "_", extra = "merge", fill = "left") %>%
  subset(., measure %in% c(2,3)) %>%
  select(., c('task', 'PC', 'term', 'estimate')) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  `colnames<-`(c('task', 'PC', 'sigma2_b', 'sigma2_w')) %>%
  mutate(label = gsub('_', '-', task)) %>%
  mutate(component = case_when(
    PC == 'PC1' ~ 'Episodic Knowledge',
    PC == 'PC2' ~ 'Intrusive Distraction',
    PC == 'PC3' ~ 'Deliberate Task-Focus',
    PC == 'PC4' ~ 'Sensory Engagement'
    ))
  
options(warn=0)

pertask_fieldmap <- rex_plot.var.field(taskfield_colored,
                 plot.point = F, plot.density = F, show.icc_slope = F, alpha.density = 0.3,
                 axis.min = 0, axis.max = 1.75) +
  geom_point(aes(x = sigma2_w, y = sigma2_b,
                 fill = color), color = "black", pch = 21, size = 4, show.legend = F) +
  facet_wrap(~ component, ncol = 2, nrow = 2) +
  ylab('B/w-Subject Variation') +
  # facet_wrap(~ component, scales = "fixed", ncol = 4, nrow = 1) +  # Separate plots per model
  scale_fill_identity() +
  geom_abline(intercept = 0, slope = 1/3, linetype = 'dashed') +
  annotate('text', label = "0.25", x = 1.25, y = 0.6, fontface = 'bold', size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  annotate('text', label = "ICC = 0.5", x = 1.25, y = 1.6, fontface = 'bold', size = 3) +
  geom_abline(intercept = 0, slope = 3, linetype = 'dashed') +
  annotate('text', label = "0.75", x = 0.3, y = 1.6, fontface = 'bold', size = 3) +
  theme(strip.text = element_text(size = 12, face = "bold", color = "black"),
        strip.placement = "outside",
        strip.background = element_rect(fill=NA,colour=NA),
        panel.spacing=unit(0.2,"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 9, color = 'black'),
        axis.text.y = element_text(face = "bold", size = 9, color = 'black'),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_blank())


```

### Example Field Map

```{r warning=FALSE}

# build a grid of variance values
grid <- expand_grid(
  sigma_w = seq(0, 1, length.out = 400),
  sigma_b = seq(0, 1, length.out = 400)
) %>%
  mutate(
    ICC = if_else(sigma_b + sigma_w == 0, NA_real_,
                   sigma_b / (sigma_b + sigma_w))
  )

# example points
pts <- tibble(
  label = c("Task-driven\n(low between,\nlow within)",
            "Low stability\n(high within,\nlow between)",
            "Subject-driven\n(high between,\nlow within)"),
  sigma_w = c(0.05, 0.45, 0.1),
  sigma_b = c(0.10, 0.10, 0.55)
)

# segment lines for each ICC-band
iccline <- function(val) {
  lab <- paste("ICC=", as.character(icc_labs[val]), sep="")
  
  if (val <= 5) {
    return(annotate(geom="text", x=0.9, y=icc_slope[val]*0.9,
                    label=lab, color="black", fontface = "bold",
                    angle=icc_angle[val], size=3.5)) 
  }
  
  return(annotate(geom="text", x=0.9/icc_slope[val], y=0.9,
                           label=lab, color="black", fontface = "bold",
                           angle=icc_angle[val], size=3.5))
}

icc_labs=seq(.1, .9, by=0.10)

icc_slope=array(0, 9); i=1;
for (z in seq(.1,.9,by=0.1)){icc_slope[i]=z/(1-z); i=i+1}

get.angle <- function(slope) {return((atan(slope))*180/pi)}

icc_angle=array(0, 9); i=1;
for (slope in icc_slope){icc_angle[i]=get.angle(slope); i=i+1}


fldmap <- ggplot(grid, aes(sigma_w, sigma_b, z = ICC)) +
  # greyscale ICC bands
  geom_contour_filled(
    breaks = seq(0, 1, by = 0.1),
    alpha = 0.7,
    show.legend = F
  ) +
  scale_fill_grey(start = 0.9, end = 0.2, name = "ICC") +
  # points
  geom_point(data = pts, aes(sigma_w, sigma_b),
    inherit.aes = FALSE, size = 3,
    shape = 24, fill = "black"
  ) +
  # descriptive labels
  geom_text(
    data = pts,
    aes(sigma_w, sigma_b, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = -0.5,
    size = 3,
    fontface = "bold"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    x = "Within-individual variation (σ²w)",
    y = "Between-individual variation (σ²b)"
  ) +
  theme_sjplot(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title.x = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )

for (i in seq(1:9)) {
  fldmap<-fldmap+iccline(i)
}

fldmap
```

## Bootstrapped pc-score ~ ICC coefficients

### Generate ICC x PC-score matrix

```{r}

genICC <- function(dat) {
  
  fm <- lme4::lmer(score ~ 1 + (1|subID),
    data = dat,
    REML = TRUE
  )
  
  output <- summary(fm)
  sigma2_r = as.numeric(output$sigma^2)
  sigma2_b = as.numeric(output$varcor$subID)
  icc1R = sigma2_b/(sigma2_b + sigma2_r)
  
  return(icc1R)
}
  
gen_table <- function(resamp) {
  
  ICC.table <- resamp %>%
    group_by(commspsych, component) %>%
    group_modify(~ tibble(
      ICC   = genICC(.x),
      score = mean(.x$score, na.rm = TRUE)
      )) %>%
    ungroup() %>%
    group_by(component, .drop = F) %>%
    mutate(ICC=atanh(ICC)) %>%
    ungroup()
  
  return(ICC.table)

}
  
```


### Resampling relationships

```{r}

boot_components <- function(dat, B=1000) {
  parts <- unique(dat$subID)
  comps <- unique(dat$component)

  out <- vector("list", B)

  for (b in seq_len(B)) {
    # resample participants with replacement
    samp_parts <- sample(parts, length(parts), replace = TRUE)

    d_b <- map_dfr(samp_parts, ~ dat[dat$subID == .x, ])

    # recompute task × component summaries
    tc_b <- gen_table(d_b)

    # fit separate model for each component
    res_b <- map_dfr(comps, function(cmp) {
      d_cmp <- tc_b %>% filter(component == cmp)

      fit <- lm(ICC ~ score, data = d_cmp)

      tibble(
        iter = b,
        Component = cmp,
        slope = coef(fit)[["score"]],
        intercept = coef(fit)[["(Intercept)"]]
      )
    })

    out[[b]] <- res_b
  }

  bind_rows(out)
}

pc_long <- pc_scores %>%
  select(subID, commspsych, PC1, PC2, PC3, PC4) %>%
  pivot_longer(PC1:PC4, names_to = "component", values_to = "score")

set.seed(123)
full_boot <- boot_components(pc_long, B = 1000)

write.csv(full_boot, 'resampled_iccPC.csv', row.names = F)

```


### Summary table

```{r}

summary_tbl <- full_boot %>%
  group_by(Component) %>%
  summarise(
    median_slope = median(slope),
    lo = quantile(slope, .025),
    hi = quantile(slope, .975),
    prop_positive = mean(slope > 0),
    prop_negative = mean(slope < 0),
    .groups = "drop"
  )

summary_tbl

```

### Plotting

```{r message=FALSE}

## wrangle plotting data
pc_long <- master_set %>%
  select(subID, commspsych, PC1, PC2, PC3, PC4) %>%
  pivot_longer(PC1:PC4, names_to = "component", values_to = "score")

pc_avs <- gen_table(pc_long)

taskcols <- task_avs_colored %>%
  select(., commspsych, task_cat, color, category_fill)

colnames(pc_avs)[2] <- "Component"

pc_avs <- pc_avs %>% 
  mutate(task_cat = case_when(
    str_detect(commspsych, "^2B|Math$|1B") ~ "ExecutiveWM",
    str_detect(commspsych, "0B|FingerTap|GoNoGo") ~ "Visuomotor",
    str_detect(commspsych, "Read|Documentary|SciFi") ~ "Naturalistic",
    str_detect(commspsych, "Memory|Friend|You") ~ "InternallyOriented"
  )) %>%
  mutate(Component = factor(Component,
                          levels = c("PC1", "PC2", "PC3", "PC4"),
                          labels = c("Episodic Knowledge",
                                     "Intrusive Distraction",
                                     "Deliberate Task-Focus",
                                     "Sensory Engagement"))) %>%
  merge(taskcols, by = c('task_cat', 'commspsych'))


## build an x-grid per component based on observed data
x_grid <- pc_avs %>%
  group_by(Component) %>%
  summarise(
    xmin = min(score),
    xmax = max(score),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(x = list(seq(xmin, xmax, length.out = 50))) %>%
  unnest(x)

median_fits <- full_boot %>%
  mutate(Component = factor(Component,
                            levels = c("PC1", "PC2", "PC3", "PC4"),
                            labels = c("Episodic Knowledge",
                                     "Intrusive Distraction",
                                     "Deliberate Task-Focus",
                                     "Sensory Engagement"))) %>%
  group_by(Component) %>%
  summarise(slope = median(slope),
            intercept = median(intercept)) %>%
  inner_join(x_grid, by = "Component") %>%
  mutate(yhat = intercept + slope * x)
  

# turn slopes/intercepts into predicted lines
boot_fits <- full_boot %>%
    mutate(Component = factor(Component,
                          levels = c("PC1", "PC2", "PC3", "PC4"),
                          labels = c("Episodic Knowledge",
                                     "Intrusive Distraction",
                                     "Deliberate Task-Focus",
                                     "Sensory Engagement")))%>%
  inner_join(x_grid, by = "Component") %>%
  mutate(yhat = intercept + slope * x)

ggplot(median_fits, aes(x = x, y = yhat)) +
  # median regression line
  geom_line(alpha = 1, colour = "black", linewidth = 1) +
  # observed task-level points
  geom_point(data = pc_avs,
             aes(x = score, y = ICC,
                 color = color, fill = task_cat),
             size = 3, show.legend = F) +
  # resampled regression lines
  geom_line(
    data = boot_fits,
    aes(x = x, y = yhat, group = iter),
    alpha = 0.02,
    colour = "grey40"
  ) +
  labs(
    x = "Mean component score (task-level)",
    y = expression(bold('ICC'[italic('Subject')]~' per Task'~" (Fisher "~italic(z)~"-transform)")),
  ) +
  facet_wrap(~ Component, scales = "free_x") +
  scale_color_identity() + 
  theme_sjplot() +
  theme(strip.text = element_text(size = 12, face = "bold", color = "black"),
        strip.placement = "outside",
        strip.background = element_rect(fill=NA,colour=NA),
        panel.spacing = unit(1, "lines"),
        axis.title.x = element_text(face = "bold", size = 12, color = "black"),
        axis.title.y = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(face = "bold", color = "black",
                                   size = 10),
        axis.text.y = element_text(face = "bold", color = "black",
                                   size = 10))


```

## Associating task-level ICCs with gradients

### Generate ICC x Gradient Matrix

```{r}

grad_coords <- grad_coords %>%
  select(., c('Task_name', contains('gradient'))) %>%
  mutate(across(gradient1:gradient5, scale)) %>%
  mutate(gradient1 = gradient1[,1],
         gradient2 = gradient2[,1],
         gradient3 = gradient3[,1],
         gradient4 = gradient4[,1],
         gradient5 = gradient5[,1])

grad_table <- function(resamp, comp=NULL, grad=NULL) {
  
  ICC.table <- resamp %>%
    group_by(commspsych, component) %>%
    # subset(., component == comp) %>%
    group_modify(~ tibble(ICC = genICC(.x))) %>%
    ungroup() %>%
    group_by(component, .drop = F) %>%
    mutate(ICC=atanh(ICC)) %>%
    ungroup() %>%
    # pivot_wider(names_from=component, values_from=ICC) %>%
    merge(grad_coords, by.x='commspsych', by.y='Task_name')
    # subset(., gradient == grad)
  
  return(ICC.table)

}

```

### Resampling relationships

```{r}

boot_gradients <- function(dat, B=1000) {
  parts <- unique(dat$subID)
  comps <- unique(dat$component)

  out <- vector("list", B)

  for (b in seq_len(B)) {
    # resample participants with replacement
    samp_parts <- sample(parts, length(parts), replace = TRUE)

    d_b <- map_dfr(samp_parts, ~ dat[dat$subID == .x, ])

    # recompute task × component summaries
    tc_b <- grad_table(d_b)

    # fit separate model for each component
    res_b <- map_dfr(comps, function(cmp) {
      d_cmp <- tc_b %>% filter(component == cmp)

      fit <- lm(ICC ~ gradient1 + gradient2 + gradient3 + gradient4 + gradient5,
                data = d_cmp)

      tibble(
        iter = b,
        Component = cmp,
        g1 = coef(fit)[["gradient1"]],
        g2 = coef(fit)[["gradient2"]],
        g3 = coef(fit)[["gradient3"]],
        g4 = coef(fit)[["gradient4"]],
        g5 = coef(fit)[["gradient5"]],
        intercept = coef(fit)[["(Intercept)"]]
      )
    })

    out[[b]] <- res_b
  }

  bind_rows(out)
}

set.seed(123)

B <- 1000

full.boot.reg <- boot_gradients(pc_long, B=B)

write.csv(full.boot.reg, 'gradicc_regression_z.csv', row.names = F)


```

### Summary table

```{r}

summary.tbl.reg <- full.boot.reg %>%
  pivot_longer(g1:g5, names_to = "gradient", values_to = "slope") %>%
  group_by(Component, gradient, .drop=F) %>%
  summarise(median_slope = median(slope),
            lo = quantile(slope, .00125),
            hi = quantile(slope, .99875),
            prop_positive = mean(slope > 0),
            prop_negative = mean(slope < 0),
            .groups = "drop"
            )

write.csv(summary.tbl.reg, 'regression_summary_z.csv', row.names = F)

```

### Plotting 

```{r}

summary.tbl.reg <- read.csv('regression_summary_z.csv')

taskcols <- task_avs_colored %>%
  select(., commspsych, task_cat, color, category_fill)

## Generate observed correlations
grad_avs <- grad_table(pc_long) %>%
  mutate(task_cat = case_when(
    str_detect(commspsych, "^2B|Math$|1B") ~ "ExecutiveWM",
    str_detect(commspsych, "0B|FingerTap|GoNoGo") ~ "Visuomotor",
    str_detect(commspsych, "Read|Documentary|SciFi") ~ "Naturalistic",
    str_detect(commspsych, "Memory|Friend|You") ~ "InternallyOriented"
  )) %>%
  merge(taskcols, by = c('task_cat', 'commspsych'))

colnames(grad_avs)[3] = "Component"

# build an x-grid per component based on observed data
grad_long <- grad_avs %>%
  pivot_longer(gradient1:gradient5,
               names_to='grad',
               values_to='coordinate') %>%
  mutate(grad = factor(grad,
                           levels = c('gradient1',
                                      'gradient2',
                                      'gradient3',
                                      'gradient4',
                                      'gradient5'),
                           labels = c('Uni. – Hetero.',
                                      'Mot. – Vis.',
                                      'DMN – Cont.',
                                      'DAN – S.Att.',
                                      'S.Mot.–V.Sem.')))

sum.tbl.plt <- summary.tbl.reg %>%
  mutate(Component = factor(Component,
                       levels = c("PC1", "PC2", "PC3", "PC4"),
                       labels = c("Episodic Knowledge",
                                     "Intrusive Distraction",
                                     "Deliberate Task-Focus",
                                     "Sensory Engagement")))

gradCOLS <- c('#ff5a5f', '#ffb400', '#007a87', '#8ce071', '#7b0051')

ggplot(sum.tbl.plt, aes(x=gradient, y=median_slope,
                            ymin=lo, ymax=hi, color = gradient)) +
  geom_pointrange(size = 0.8, linewidth = 1) +
  geom_hline(yintercept=0, lty=2, color='black', size = 1) +
  scale_x_discrete(name="",
                   limits = c("g5", "g4", "g3", "g2", "g1"),
                   labels = rev(c('Uni. – Hetero.',
                                      'Mot. – Vis.',
                                      'DMN – Cont.',
                                      'DAN – S.Att.',
                                      'S.Mot.–V.Sem.'))) +
  scale_y_continuous(name="Median slope (95% CI)") +
  scale_color_manual(name = "Gradient",
                     values = gradCOLS,
                     limits = rev(c("g5", "g4", "g3", "g2", "g1")),
                     labels = c('Uni. – Hetero.',
                                'Mot. – Vis.',
                                'DMN – Cont.',
                                'DAN – S.Att.',
                                'S.Mot.–Lim.VAN')) +
  coord_flip() +
  facet_wrap(~ Component, scales = "fixed", ncol=4, nrow=1) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        # axis.text.y = element_text(face = "bold.italic", size = 10, color = "black"),
        axis.title.x = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(face = "bold", size = 10, color = "black"),
        strip.text = element_text(size = 12, face = "bold", color = "black"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.title = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(fill=NA,colour=NA),
        panel.spacing=unit(0.3,"cm"))

```

### MDN Similarity

```{r}
mdnsim <- as.data.frame(t(read.csv('regression_mdnsimilarity_z.csv'))) %>%
  `colnames<-`(c('PC1', 'PC2', 'PC3', 'PC4')) %>%
  pivot_longer(PC1:PC4,
               names_to='comp',
               values_to='corr') %>%
  group_by(comp) %>%
  summarise(
    mean = mean(corr),
    median = median(corr),
    lo = quantile(corr, .025),
    hi = quantile(corr, .975),
    prop_positive = mean(corr > 0),
    prop_negative = mean(corr < 0),
    .groups = "drop"
  ) %>%
  mutate(label = c(
        "Episodic Knowledge",
        "Intrusive Distraction",
        "Deliberate Task-Focus",
        "Sensory Engagement"
      )) %>%
  mutate(label = factor(label, levels = c(
        "Episodic Knowledge",
        "Intrusive Distraction",
        "Deliberate Task-Focus",
        "Sensory Engagement"
      )
    )
  ) %>%
  droplevels()

dotCOLS = c("#a6d8f0","#f9b282","#a7d876","#f5f496")
barCOLS = c("#008fd5","#de6b35","#65ac1f", "#f3f00f")

ggplot(mdnsim,
       aes(x=label, y=mean, ymin=lo, ymax=hi, color = label, fill = label)) +
  geom_linerange(size=5, show.legend = F) +
  geom_hline(yintercept=0, lty=2, color='black', size = 1) +
  scale_fill_manual(values=barCOLS) +
  scale_color_manual(values=dotCOLS) +
#specify position here too
  geom_point(size=4, shape=21, color = "white", stroke = 0.5, show.legend = F) +
  scale_x_discrete(name="Stability per Thought-Pattern",
                   limits=rev(c("Episodic Knowledge",
                            "Intrusive Distraction",
                            "Deliberate Task-Focus",
                            "Sensory Engagement"))) +
  scale_y_continuous(name="Similarity to MDN (Bootstrapped 95% CI)") +
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold.italic", size = 10, color = "black"),
        axis.title.x = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(face = "bold", size = 10, color = "black"))

```


