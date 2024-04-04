
# Script to reproduce the results and figures from the paper: Longitudinal dynamics between the central nodes in the symptoms network of Borderline Personality Disorder: an intraindividual network analysis.

# April 2024 
# by: Said Jiménez (said.ejp@gmail.com)

# Used packages
pacman::p_load(tidyverse, bootnet, qgraph, psychonetrics)

# Data can be found at: https://dataverse.nl/dataset.xhtml?persistentId=doi:10.34894/NVROJB 
load("dataverse_files/vonKlipstein2021.RData")

# Preprocessing -------
data_bpd <- 
  tibble(shareD) %>% 
  filter(!is.na(BPDSIsum))

times <- c(0, 6, 12, 18, 24)

indices <-
  data_bpd %>% 
  filter(study %in% c("Giesen-Bloo", "GroupST"), 
         time %in% times) %>% 
  count(ID) %>% 
  filter(n == length(times)) %>% 
  pull(ID)

indices 

data_bpd <- 
  data_bpd %>% 
  filter(ID %in% indices, 
         time %in% times) 

data_bpd %>% 
  distinct(ID)


data_bpd %>% 
  count(time)

data_bpd %>% 
  count(study)

data_bpd %>% 
  count(study, time)

# Selected nodes -----
# read our previous publication: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0289101
items_av <- 
  data_bpd %>% 
  select(ID, contains(c("6.5", "3.2", "4.8", "5.9", "7.1", "8.3", "9.8", "1.6", "2.1")))

vars <- colnames(items_av)[-1]

items_av <- 
  items_av %>% 
  select(ID, vars[order(vars)])


vars_ord <- colnames(items_av)[-1]

# Fitting model -----
model <- 
  ml_gvar(
    data = data.frame(items_av), 
    vars = vars_ord, 
    idvar = "ID", 
    standardize = "z"
  ) %>% 
  runmodel


model2 <- 
  model %>% 
  prune(
    alpha = 0.01, 
    recursive = FALSE
  ) %>% 
  modelsearch(verbose = FALSE)

# Fit measures:
model2 %>% fit

compare(model, model2)

# Extract networks:
labels <- str_remove(colnames(items_av)[-1], pattern = "BPDSI")

nets_dets <- 
  map(
    .x = c("PDC", "omega_zeta_within", "omega_zeta_between"), 
    .f = ~ getmatrix(model2, .x)
  )

nodes <- c("Feeling angry",
           "Identity change",
           "Binge eating", 
           "Want to commit suicide",
           "Feeling empty",
           "Scream and slam doors",
           "Unfair treatment", 
           "Hear love from someone", 
           "Unstable relationships")

nodes_ord <- nodes[order(vars)]

lay <- averageLayout(nets_dets[[1]],  layout = "spring")
qgraph(nets_dets[[1]], layout = lay, labels = labels, 
              title = "Temporal", shape = "circle",

       
              #vsize = 10, vsize2 = 8,
              mar = rep(8, 4), 
              #asize = 7, 
              edge.labels = F,
              theme = "colorblind", 
              nodeNames = nodes_ord,
              legend.cex = .5,
              label.color = 'black'
              #filetype = "png", 
              #filename = "temporal"
       )

# Figures 1 and 2 in main manuscript
nets_dets_graph <-
  map2(
    .x = nets_dets, 
    .y = c("A. Temporal", "B. Contemporaneous", "C. Between-persons"),
    .f = function(i, j) {
      
      png(paste0(j, ".png"), units = "cm", height = 10, width = 15, res = 600)
      
      qgraph(i, layout = lay, labels = labels, 
             shape = "circle", 
             title = j,
             #vsize = 10, vsize2 = 8,
             mar = rep(8, 4), 
             #asize = 7, 
             edge.labels = F,
             theme = "colorblind",
             nodeNames = nodes_ord, legend.cex = .35
             #filetype = "png", 
             #filename = j
      )
      
      dev.off()
      
    }
  )


# Bootstrapping -----------------------------------------------------------



set.seed(2024)

n_subjects <- 
  items_av %>% 
  select(ID) %>% 
  distinct() %>% 
  nrow()

n_reps <- 1000

out <- list()

for (i in 1:n_reps) {
  
  out[[i]] <- 
    items_av %>% 
    select(ID) %>% 
    distinct() %>% 
    sample_n(
      size = n_subjects, 
      replace = TRUE
    ) %>% 
    pull(ID)
  
}


test <- 
  map(
  .x = 1:n_reps, 
  .f = function(i) { 
    
    id_count <- 
      tibble(
        ID = out[[i]]
      ) %>% 
      count(ID) 
    
    id_nest <-
      items_av %>% 
      group_by(ID) %>% 
      nest() %>% 
      filter(ID %in% id_count$ID) 
    
    df_boot <- 
      id_nest %>% 
      inner_join(id_count) %>% 
      uncount(n, .id = "reps") %>% 
      unnest(data) %>% 
      ungroup() %>% 
      mutate(ID_2 = rep(1:n_subjects, each = 5))
    
    print(paste("Rep", i, "started!"))
    
    
    model_boot <- 
      ml_gvar(
        data = data.frame(df_boot), 
        vars = vars, 
        idvar = "ID_2", 
        standardize = "z"
      ) %>% 
      runmodel
    
    }
)

# write_rds(test, file = "boot_1000_replicates_var_bpd_2.rds", compress = "gz")
  
# test <- read_rds("boot_1000_replicates_var_bpd_2.rds")

boot_temporal_det <- map(test, ~ getmatrix(.x, "PDC"))
boot_contemp_det <- map(test, ~ getmatrix(.x, "omega_zeta_within"))
boot_between_det <- map(test, ~ getmatrix(.x, "omega_zeta_between"))


node_names_short <- str_remove(vars, pattern = "BPDSI")
  
boot_det_tidy_df <- 
  map(
    .x = list(boot_temporal_det, boot_contemp_det, boot_between_det), 
    .f = function(lista) { 
      
      map_df(
        .x = 1:n_reps, 
        .f = function(i) { 
          
          cor_ex <- lista[[i]]
          
          colnames(cor_ex) <- node_names_short
          
          cor_ex %>% 
            as_tibble() %>% 
            mutate(term_a = node_names_short) %>% 
            select(term_a, everything()) %>% 
            pivot_longer(
              cols = 2:last_col(), 
              names_to = "term_b", 
              values_to = "coef"
            )
        },
        .id = "replicate"
      )
    }
  )



names(boot_det_tidy_df) <- c("Temporal", "Contemporaneous", "Between_person")    


fig_boot_temp_det <- 
  boot_det_tidy_df$Temporal %>% 
  unite("A_B", term_a:term_b, sep = " -> ") %>% 
  group_by(A_B) %>% 
  summarize(
    mean = mean(coef, na.rm = T), 
    low = quantile(coef, probs = .025, na.rm = T), 
    high = quantile(coef, probs = .975, na.rm = T)
  ) %>% 
  mutate(A_B = fct_reorder(A_B, mean, .fun = max)) %>% 
  top_n(30) %>% 
  ggplot(aes(x = A_B, y = mean)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
  geom_pointrange(aes(ymin = low, ymax = high), 
                  col = "#203DD9", size = 1/8) +
  scale_y_continuous(limits = c(-1, 1)) +
  #geom_point(
  #  data = temporal_full,
  #  aes(x = A_B, y = mean_orig), 
  #  col = "red"
  #  ) +
  labs(x = NULL, 
       y = "Estimate", 
       tag = "A") +
  coord_flip() +
  theme_light() +
  facet_wrap(~"Temporal") +
  theme(strip.text = element_text(face = "bold", color = "black"), 
        axis.text.y = element_text(size = 8))

fig_boot_temp_det

fig_boot_cont_det <-
  boot_det_tidy_df$Contemporaneous %>% 
  unite("A_B", term_a:term_b, sep = " - ") %>% 
  group_by(A_B) %>% 
  summarize(
    mean = mean(coef, na.rm = T), 
    low = quantile(coef, probs = .025, na.rm = T), 
    high = quantile(coef, probs = .975, na.rm = T),
  ) %>% 
  filter(mean != 0) %>% 
  distinct(mean, .keep_all = T) %>% 
  mutate(A_B = fct_reorder(A_B, mean, .fun = max)) %>% 
  ggplot(aes(x = A_B, y = mean)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
  geom_pointrange(aes(ymin = low, ymax = high), 
                  col = "#203DD9", size = 1/8) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = NULL, 
       y = "Estimate", 
       tag = "B") +
  coord_flip() +
  theme_light() +
  facet_wrap(~"Contemporaneous") +
  theme(strip.text = element_text(face = "bold", color = "black"),
        axis.text.y = element_text(size = 7))

fig_boot_cont_det


fig_boot_between_det <-
  boot_det_tidy_df$Between_person %>% 
  unite("A_B", term_a:term_b, sep = " - ") %>% 
  group_by(A_B) %>% 
  summarize(
    mean = mean(coef, na.rm = T), 
    low = quantile(coef, probs = .025, na.rm = T), 
    high = quantile(coef, probs = .975, na.rm = T),
  ) %>% 
  filter(mean != 0) %>% 
  distinct(mean, .keep_all = T) %>% 
  mutate(A_B = fct_reorder(A_B, mean, .fun = max)) %>% 
  ggplot(aes(x = A_B, y = mean)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
  geom_pointrange(aes(ymin = low, ymax = high), 
                  col = "#203DD9") +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = NULL, 
       y = "Estimate", 
       tag = "C") +
  coord_flip() +
  theme_light() +
  facet_wrap(~"Between-persons") +
  theme(strip.text = element_text(face = "bold", color = "black"))

fig_boot_between_det

fig_boot <-
  gridExtra::grid.arrange(fig_boot_temp_det, 
                          fig_boot_cont_det,
                          ncol = 2)



boot_det_tidy_df$Temporal %>% 
  unite("A_B", term_a:term_b, sep = " -> ") %>% 
  group_by(A_B) %>% 
  summarize(
    mean = mean(coef, na.rm = T), 
    low = quantile(coef, probs = .025, na.rm = T), 
    high = quantile(coef, probs = .975, na.rm = T),
  ) %>% filter(low > 0) %>% arrange(desc(mean))

# Figure 4 in main manuscript
ggsave(
  filename = "fig_stability.png", 
  path = "figures_temporal_bpd/",
  plot = fig_boot, 
  device = "png", 
  dpi = 350, 
  units = "cm", 
  width = 15, 
  height = 10
)

# Centrality nodes ----------------------------------------------------

names(nets_dets) <- c("Temporal", "Contemporaneous", "Between_person")

df_cent <- 
  map(
  .x = nets_dets, 
  .f = ~ centralityTable(.x, standardized = TRUE)
) %>% 
  bind_rows(.id = "network") %>% 
  select(-graph, -type) %>% 
  mutate(node = factor(node, levels = 1:9, labels = labels)) %>% 
  as_tibble()

fig_cent <- 
df_cent %>% 
  filter(network != "Between_person", 
         str_detect(measure, "Expected", negate = T)) %>% 
  split(.$network) %>% 
  map(
    .f = function(df) {
      
      df %>% 
        ggplot(aes(x = node, y = value)) +
        geom_hline(yintercept = 0, lty = 2, alpha = .5) +
        geom_segment(aes(xend = node, y = 0, yend = value), col = "#203DD9") +
        geom_point(col = "#203DD9") +
        scale_y_continuous(limits = c(-3, 3)) +
        labs(x = "Node", 
             y = "Z-score", 
             title = df[1, 1]) +
        facet_grid(~ measure) +
        coord_flip() +
        theme_light() +
        theme(strip.text = element_text(color = "black"), 
              plot.title = element_text(hjust = .5, face = "bold"))

      
    }
  )

figs_cent <- 
  gridExtra::grid.arrange(
  fig_cent[[2]] + labs(tag = "A"), 
  fig_cent[[1]] + labs(tag = "B")
)

# Figure 3 in main manuscript
ggsave(
      plot = figs_cent,
      filename = "fig_centrality.png", 
      path = "figures_temporal_bpd/", 
      device = "png", 
      dpi = 350, 
      units = "cm", 
      width = 17, 
      height = 15
      )


#map2(
#  .x = fig_cent, 
#  .y = c("fig_cent_contemporaneous.png", "fig_cent_temporal.png"),
#  .f = function(p, n) {
#    
#    
#    ggsave(
#      plot = p,
#      filename = n, 
#      path = "figures_temporal_bpd/", 
#      device = "png", 
#      dpi = 350, 
#      units = "cm", 
#      width = 15, 
#      height = 10
#      )
#    
#    
#  }
#)


# Standard errors ---------------------------------------------------------

# Table 1 in main manuscript
matrix_coefs <- 
  boot_det_tidy_df %>% 
  map(
    .f = function(df) {
      
      df_res <- 
      df %>% 
        group_by(term_a, term_b) %>% 
        summarize(
          r = round(mean(coef), 2),
          se = round(sd(coef), 2), 
          .groups = "drop"
        ) %>% 
        mutate(
          r = ifelse(r == 0, NA, r), 
          r = str_replace(r, pattern = "^0", ""),
          r =  str_replace(r, pattern = "^-0", "-"),
          se = ifelse(se == 0, NA, se),
          se = str_replace(se, pattern = "^0", ""),
          r_se = paste0(r, " [", se, "]"), 
          r_se = ifelse(str_detect(r_se, pattern = "NA"),  NA, r_se)
          )
       
      
      df_res %>% 
        pivot_wider(
          id_cols = term_a, 
          names_from = term_b, 
          values_from = r_se
        )

    }
  )

