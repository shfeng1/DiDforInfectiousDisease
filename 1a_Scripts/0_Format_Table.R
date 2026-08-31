# Summarize power and bias results. Print into a nice table
format.tbl <- function(power.df, bias.original.df, bias.AME.df) {
  typeI.error <- power.df %>% filter(eff.multi==1) %>% group_by(model) %>%
    summarise(p = round(mean(p), 1)) %>% dplyr::select(model, p) %>%
    arrange(match(model, model.list))
  power <- power.df %>% filter(eff.multi==1.1) %>% group_by(model) %>%
    summarise(power = round(mean(p), 1)) %>% dplyr::select(model, power) %>%
    arrange(match(model, model.list))
  bias.original <- bias.original.df %>% group_by(model) %>% 
    summarise(mean = format(round(mean(abs(eff.bias)), 2), nsmall=2), 
              min = format(round(min(eff.bias), 2), nsmall=2), 
              max = format(round(max(eff.bias), 2), nsmall=2)) %>%
    mutate(bias_original = paste0(mean, " (", min, ", ", max, ")")) %>%
    arrange(match(model, model.list))
  bias.AME <- bias.AME.df %>% group_by(model) %>% 
    filter(model != "true") %>%
    summarise(mean = format(round(mean(abs(bias.fit)), 1), nsmall=1),
              min = format(round(min(bias.fit), 1), nsmall=1),
              max = format(round(max(bias.fit), 1), nsmall=1)) %>%
    mutate(bias_AME = paste0(mean, " (", min, ", ", max, ")")) %>%
    arrange(match(model, model.list))
  
  typeI.error %>% cbind(power$power) %>% cbind(bias.original$bias_original) %>% cbind(bias.AME$bias_AME)
}

format.tbl.A3 <- function(bias.original.df, bias.AME.df) {
  bias.original <- bias.original.df %>%
    group_by(model) %>%
    summarise(mean = format(round(mean(abs(eff.bias)), 2), nsmall=2),
              min = format(round(min(eff.bias), 2), nsmall=2),
              max = format(round(max(eff.bias), 2), nsmall=2)) %>%
    mutate(bias_original = paste0(mean, " (", min, ", ", max, ")")) %>%
    dplyr::select(model, bias_original)
  bias.AME <- bias.AME.df %>%
    group_by(model) %>%
    summarise(mean = format(round(mean(abs(bias.fit)), 1), nsmall=1),
              min = format(round(min(bias.fit), 1), nsmall=1),
              max = format(round(max(bias.fit), 1), nsmall=1)) %>%
    mutate(bias_AME = paste0(mean, " (", min, ", ", max, ")")) %>%
    dplyr::select(model, bias_AME)
  bias.AME.adj <- bias.AME.df %>%
    group_by(model) %>%
    summarise(mean = format(round(mean(abs(bias.adj)), 1), nsmall=1),
              min = format(round(min(bias.adj), 1), nsmall=1),
              max = format(round(max(bias.adj), 1), nsmall=1)) %>%
    mutate(bias_AME_correct = paste0(mean, " (", min, ", ", max, ")")) %>%
    dplyr::select(model, bias_AME_correct)

  bias.original %>%
    left_join(bias.AME, by="model") %>%
    left_join(bias.AME.adj, by="model")
}
