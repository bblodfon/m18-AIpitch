# M18 PANCAIM Meeting
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(rpart.plot)
set.seed(42)

# Task Lung ----
?survival::lung

# Data Preprocessing
lung = survival::lung
lung %>% as_tibble()

lung$status = (lung$status == 2L) # 2 is death so convert to 1
lung = lung %>% select(-inst) # remove Institution code (irrelevant for us)
lung$ph.ecog = as.integer(lung$ph.ecog)

task = as_task_surv(x = lung, time = 'time', event = 'status', id = 'lung')
task$missings() # missing values!

# Impute NA's
impute_po = po("imputelearner", lrn('regr.rpart'))
task_lung = impute_po$train(list(task))[[1]]
task_lung$missings() # no missing values

# split dataset to train and test
train_set = sample(task_lung$nrow, 0.8 * task_lung$nrow)
test_set = setdiff(seq_len(task_lung$nrow), train_set)
intersect(train_set, test_set)

# Learners ----
cox = lrn('surv.coxph', id = 'CoxPH')
rpart = lrn('surv.rpart', id = 'SurvivalTree')
ranger = lrn('surv.ranger', verbose = FALSE, id = 'SurvivalForest')
#ctree = lrn('surv.ctree')

learners = list(cox, rpart, ranger)

# Train & Test ----
preds  = list()
cindex = list()
for (learner in learners) {
  pred = learner$train(task_lung, row_ids = train_set)$
    predict(task_lung, row_ids = test_set)
  preds[[learner$id]]  = pred
  cindex[[learner$id]] = pred$score()
}

# Benchmark barplot ----
data_tbl = bind_rows(cindex, .id = 'learner')

data_tbl %>%
  ggplot(aes(x = learner, y = surv.cindex, fill = learner)) +
  geom_bar(stat = "identity") + # show.legend = FALSE
  scale_fill_brewer(palette = 'Set1') +
  ggpubr::theme_classic2() +
  ylab('C-index') +
  coord_flip() +
  geom_hline(yintercept = 0.5, linetype = 'dotted') +
  geom_text(x = 2.5, y = 0.47, label = "Random", angle = 90) +
  theme(axis.title.y = element_blank())
ggsave(filename = 'img/cmp_learners.png', width = 6, height = 2)

# Plot tree ---
fit = rpart$model

png(filename = 'img/tree.png', units = 'in', width = 6, height = 6, res = 300)
rpart.plot(fit, type = 2, extra = 101, box.palette = "GnBu", branch.lty = 3,
  shadow.col = "gray", nn = FALSE, digits = 3)
dev.off()

# Plot Importance scores ----
var_imp = tibble::enframe(rpart$importance(), name = 'Variable', value = 'Importance')
var_imp

p = var_imp %>% ggpubr::ggdotchart(x = 'Variable', y = 'Importance',
  add = 'segment', dot.size = 5, rotate = TRUE, color = 'Variable',
  xlab = 'Features', sorting = 'ascending')
ggpubr::ggpar(p, legend = 'none')
ggsave(filename = 'tree_varimp.png', width = 4, height = 4)

# Compose Survival Distributions ----
kaplan_pred = lrn('surv.kaplan')$train(task_lung, row_ids = train_set)$
  predict(task_lung, row_ids = test_set)
nelson_pred = lrn('surv.nelson')$train(task_lung, row_ids = train_set)$
  predict(task_lung, row_ids = test_set)

# check: different survival distr
stopifnot(kaplan_pred$distr$survival(500) != nelson_pred$distr$survival(500))

# distr compose pipe operators
po_distrcompose_aft = po("distrcompose", param_vals = list(form = 'aft', overwrite = TRUE))
po_distrcompose_ph  = po("distrcompose", param_vals = list(form = 'ph', overwrite = TRUE))

# Compose distr survival tree predictions of type `aft`, using Kaplan-Meier base predictions
preds_km_aft     = po_distrcompose_aft$predict(list(base = kaplan_pred, pred = preds$SurvivalTree))[[1]]
# Compose distr survival tree predictions of type `aft`,  using Nelson-Aalen base predictions
preds_nelson_aft = po_distrcompose_aft$predict(list(base = nelson_pred, pred = preds$SurvivalTree))[[1]]

# Compose distr survival tree predictions of type `ph`, using Kaplan-Meier base predictions
preds_km_ph     = po_distrcompose_ph$predict(list(base = kaplan_pred, pred = preds$SurvivalTree))[[1]]
# Compose distr survival tree predictions of type `ph`,  using Nelson-Aalen base predictions
preds_nelson_ph = po_distrcompose_ph$predict(list(base = nelson_pred, pred = preds$SurvivalTree))[[1]]

# check: crank didn't change
stopifnot(preds_km_aft$crank == preds_nelson_aft$crank,
          preds_km_aft$crank == preds_km_ph$crank,
          preds_km_aft$crank == preds_nelson_ph$crank)

# check: survival values are all different (at least in a specific time point)
stopifnot(preds_km_aft$distr$survival(c(500)) != preds_nelson_aft$distr$survival(c(500)))
stopifnot(preds_km_aft$distr$survival(c(500)) != preds_km_ph$distr$survival(c(500)))
stopifnot(preds_km_aft$distr$survival(c(500)) != preds_nelson_ph$distr$survival(c(500)))
stopifnot(preds_nelson_aft$distr$survival(c(500)) != preds_km_ph$distr$survival(c(500)))
stopifnot(preds_nelson_aft$distr$survival(c(500)) != preds_nelson_ph$distr$survival(c(500)))
stopifnot(preds_km_ph$distr$survival(c(500)) !=  preds_nelson_ph$distr$survival(c(500)))

# make a list of the composed distr predictions
distr_preds = list(kaplan_aft = preds_km_aft, nelson_aft = preds_nelson_aft,
  kaplan_ph = preds_km_ph, nelson_ph = preds_nelson_ph)

# Test case - 'bad' individual ----
# patient's data (died early)
ind = 42
ind_data = task_lung$data(rows = test_set)[ind]
ind_data

# relative risk score
preds$SurvivalTree$crank[ind]

# find probability of survival for different time points (1, 2 and 3 years)
days_1y = 365
days_2y = 2*365
days_3y = 3*365
times = c(days_1y, days_2y, days_3y)

# base estimator is not important (kaplan or nelson) => form is what matters!!!
ind_res = bind_rows(
  preds_km_aft$distr$survival(times)[,ind] %>%
    enframe(name = 'times', value = 'surv_prob') %>%
    add_column(type = 'kaplan_aft'),
  preds_nelson_aft$distr$survival(times)[,ind] %>%
    enframe(name = 'times', value = 'surv_prob') %>%
    add_column(type = 'nelson_aft'),
  preds_km_ph$distr$survival(times)[,ind] %>%
    enframe(name = 'times', value = 'surv_prob') %>%
    add_column(type = 'kaplan_ph'),
  preds_nelson_ph$distr$survival(times)[,ind] %>%
    enframe(name = 'times', value = 'surv_prob') %>%
    add_column(type = 'nelson_ph'))
ind_res

# plot the 4 different survival curves
surv_data = sapply(distr_preds, function(pred) {
  # survival probabilities on the training set time points
  1 - pred$distr$getParameterValue('cdf')[ind,]
}) %>% as_tibble(rownames = 'times') %>% mutate(times = as.integer(times))

# pick some survival probabilities for annotation
sp_1y_ph  = ind_res %>% filter(times == days_1y, type == 'kaplan_ph') %>% pull(surv_prob)
sp_1y_aft = ind_res %>% filter(times == days_1y, type == 'kaplan_aft') %>% pull(surv_prob)
sp_2y_ph  = ind_res %>% filter(times == days_2y, type == 'kaplan_ph') %>% pull(surv_prob)
sp_2y_aft = ind_res %>% filter(times == days_2y, type == 'kaplan_aft') %>% pull(surv_prob)

surv_data %>%
  tidyr::pivot_longer(cols = starts_with('kaplan'), # keep only kaplan estimator
    names_to = 'estimator', values_to = 'surv_prob') %>%
  ggplot(aes(x = times, y = surv_prob, group = estimator)) +
  geom_line(aes(color = estimator), size = 1) +
  scale_color_brewer(palette = 'Set1', name = 'Distr Form', labels = c('AFT', 'PH')) +
  # time-of-death
  geom_vline(xintercept = ind_data$time, linetype = 'dashed', color = 'grey42') +
  annotate('text', x = ind_data$time - 50, y = 0.5, label = 'Actual time of death', angle = 90,
    size = 7) +
  # 1-year mark
  geom_vline(xintercept = days_1y, linetype = 'dashed', color = 'purple') +
  annotate('text', x = days_1y - 50, y = 0.5, label = '1 year', angle = 90,
    size = 7) +
  # AFT point
  annotate('point', x = days_1y, y = sp_1y_aft,
    colour = 'black', size = 3, shape = 4, stroke = 2) +
  annotate('text', x = days_1y + 60, y = sp_1y_aft + 0.05,
    label = round(sp_1y_aft, 2), size = 5) +
  # PH point
  annotate('point', x = days_1y, y = sp_1y_ph,
    colour = 'black', size = 3, shape = 4, stroke = 2) +
  annotate('text', x = days_1y + 70, y = sp_1y_ph + 0.05,
    label = round(sp_1y_ph, 4), size = 5) +
  # 2-year mark
  geom_vline(xintercept = days_2y, linetype = 'dashed', color = 'purple') +
  annotate('text', x = days_2y - 50, y = 0.5, label = '2 years', angle = 90,
    size = 7) +
  # AFT point
  annotate('point', x = days_2y, y = sp_2y_aft,
    colour = 'black', size = 3, shape = 4, stroke = 2) +
  annotate('text', x = days_2y + 60, y = sp_2y_aft + 0.05,
    label = round(sp_2y_aft, 2), size = 5) +
  # PH point
  annotate('point', x = days_2y, y = sp_2y_ph,
    colour = 'black', size = 3, shape = 4, stroke = 2) +
  annotate('text', x = days_2y + 70, y = sp_2y_ph + 0.05,
    label = round(sp_2y_ph, 4), size = 5) +
  labs(x = 'Days', y = 'Survival Probability') +
  mlr3viz::theme_mlr3(base_size = 14, legend = 'right')
ggsave(filename = 'img/surv_plot_bad_ind.png', width = 7, height = 5)

# Test case - 'good' individual ----
# patient's data (died at ~2 years)
ind = 14
ind_data2 = task_lung$data(rows = test_set)[ind]
ind_data2

# relative risk score
preds$SurvivalTree$crank[ind]

# find probability of survival for different time points (1, 2 and 3 years)
# base estimator is not important (kaplan or nelson) => form is what matters!!!
ind_res2 = bind_rows(
  preds_km_aft$distr$survival(times)[,ind] %>%
    enframe(name = 'times', value = 'surv_prob') %>%
    add_column(type = 'kaplan_aft'),
  preds_nelson_aft$distr$survival(times)[,ind] %>%
    enframe(name = 'times', value = 'surv_prob') %>%
    add_column(type = 'nelson_aft'),
  preds_km_ph$distr$survival(times)[,ind] %>%
    enframe(name = 'times', value = 'surv_prob') %>%
    add_column(type = 'kaplan_ph'),
  preds_nelson_ph$distr$survival(times)[,ind] %>%
    enframe(name = 'times', value = 'surv_prob') %>%
    add_column(type = 'nelson_ph'))
ind_res2

# plot the 4 different survival curves
surv_data2 = sapply(distr_preds, function(pred) {
  # survival probabilities on the training set time points
  1 - pred$distr$getParameterValue('cdf')[ind,]
}) %>% as_tibble(rownames = 'times') %>% mutate(times = as.integer(times))

# pick some survival probabilities for annotation
sp_1y_ph2  = ind_res2 %>% filter(times == days_1y, type == 'kaplan_ph') %>% pull(surv_prob)
sp_1y_aft2 = ind_res2 %>% filter(times == days_1y, type == 'kaplan_aft') %>% pull(surv_prob)
sp_2y_ph2  = ind_res2 %>% filter(times == days_2y, type == 'kaplan_ph') %>% pull(surv_prob)
sp_2y_aft2 = ind_res2 %>% filter(times == days_2y, type == 'kaplan_aft') %>% pull(surv_prob)

surv_data2 %>%
  tidyr::pivot_longer(cols = starts_with('kaplan'), # keep only kaplan estimator
    names_to = 'estimator', values_to = 'surv_prob') %>%
  ggplot(aes(x = times, y = surv_prob, group = estimator)) +
  geom_line(aes(color = estimator), size = 1) +
  scale_color_brewer(palette = 'Set1', name = 'Distr Form', labels = c('AFT', 'PH')) +
  # time-of-death
  geom_vline(xintercept = ind_data2$time, linetype = 'dashed', color = 'grey42') +
  annotate('text', x = ind_data2$time + 50, y = 0.8, label = 'Actual time of death', angle = 90,
    size = 5) +
  # 1-year mark
  geom_vline(xintercept = days_1y, linetype = 'dashed', color = 'purple') +
  annotate('text', x = days_1y - 50, y = 0.55, label = '1 year', angle = 90,
    size = 7) +
  # AFT point
  annotate('point', x = days_1y, y = sp_1y_aft2,
    colour = 'black', size = 3, shape = 4, stroke = 2) +
  annotate('text', x = days_1y + 60, y = sp_1y_aft2 + 0.05,
    label = round(sp_1y_aft2, 2), size = 5) +
  # PH point
  annotate('point', x = days_1y, y = sp_1y_ph2,
    colour = 'black', size = 3, shape = 4, stroke = 2) +
  annotate('text', x = days_1y + 70, y = sp_1y_ph2 + 0.05,
    label = round(sp_1y_ph2, 2), size = 5) +
  # 2-year mark
  geom_vline(xintercept = days_2y, linetype = 'dashed', color = 'purple') +
  annotate('text', x = days_2y - 50, y = 0.55, label = '2 years', angle = 90,
    size = 7) +
  # AFT point
  annotate('point', x = days_2y, y = sp_2y_aft2,
    colour = 'black', size = 3, shape = 4, stroke = 2) +
  annotate('text', x = days_2y + 60, y = sp_2y_aft2 + 0.05,
    label = round(sp_2y_aft2, 2), size = 5) +
  # PH point
  annotate('point', x = days_2y, y = sp_2y_ph2,
    colour = 'black', size = 3, shape = 4, stroke = 2) +
  annotate('text', x = days_2y + 70, y = sp_2y_ph2 + 0.05,
    label = round(sp_2y_ph2, 3), size = 5) +
  labs(x = 'Days', y = 'Survival Probability') +
  mlr3viz::theme_mlr3(base_size = 14, legend = 'right')
ggsave(filename = 'img/surv_plot_good_ind.png', width = 7, height = 5)
