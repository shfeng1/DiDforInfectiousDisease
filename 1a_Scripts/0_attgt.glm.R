attgt.glm <- function(yname, tname, idname, gname, weightsname, 
                      control_group = "notyettreated", data_in,
                      family) {
  y <- data_in[[yname]]
  time <- data_in[[tname]]
  ID <- as.numeric(data_in[[idname]])
  group <- data_in[[gname]]
  glist <- sort(unique(group[group>0])) # ever-treated groups
  
  if (is.null(weightsname)) { # check if weights are specified
    wt <- rep(1, length(y))
  } else {
    wt <- data_in[[weightsname]]
  }
  
  df <- data.frame(cbind(ID, group, time, y, wt))
  
  ATT_gt <- matrix(nrow = length(glist), ncol = length(min(glist):max(time)))
  rownames(ATT_gt) <- glist
  colnames(ATT_gt) <- min(glist):max(time)
  
  CI.out <- readRDS("CI.rds")
  for (g in glist) {
    for (t in g:max(time)) {
      # print(paste0("Working on g=", g, ", t=", t))
      df.gt <- df %>%
        filter( (group==g)| # treated
                  (group==0)|(group>t), # never-treated or not-yet-treated at time t
                (time==t)|(time==(g-1)) # compares periods t and (g-1)
        ) %>%
        mutate(trt_post = (group==g)*(time==t))
      fit.gt <- glm(y ~ -1+factor(ID)+factor(time)+trt_post, data=df.gt,
                    weights = df.gt$wt, family = family)
      coef <- as.numeric(tail(coef(fit.gt), 1))
      ATT_gt[as.character(g), as.character(t)] <- coef
    }
  }
    
  weight <- sapply(glist, function(g) mean(df$group[df$group>0]==g))
  out <- c()
  for (col in 1:ncol(ATT_gt)) {
    ind <- 1:length(na.omit(ATT_gt[,col]))
    out <- c(out, as.numeric(weight[ind] %*% ATT_gt[ind,col]))
  }
  
  return(out)
}