#Create a function for breakpoint estimation
find_break_rss <- function(time, freq, log = FALSE, min_seg = 20,
                           do_perm = FALSE, nperm = 1000) {
  
  if (log){
    freq <- log(freq + 1e-8)
  }
  
  n <- length(time)
  rss <- rep(NA, n)
  gain <- rep(NA, n)
  
  # global model
  fit_all <- lm(freq ~ time)
  rss_all <- sum(resid(fit_all)^2)
  
  # compute breakpoints
  for(t in min_seg:(n - min_seg)){
    
    fit1 <- lm(freq[1:t] ~ time[1:t])
    fit2 <- lm(freq[t:n] ~ time[t:n])
    
    rss[t] <- sum(resid(fit1)^2) + sum(resid(fit2)^2)
    gain[t] <- rss_all - rss[t]
  }
  
  # best breakpoint
  bp <- which.max(gain)
  
  best_gain <- gain[bp]
  rel_gain <- best_gain / rss_all
  
  # --------------------------
  # permutation test (optional)
  # --------------------------
  perm_gains <- NULL
  p_value <- NA
  
  if (do_perm) {
    perm_gains <- numeric(nperm)
    
    for (i in 1:nperm) {
      freq_perm <- sample(freq)
      
      rss_perm <- rep(NA, n)
      
      # recompute rss for permuted data
      for(t in min_seg:(n - min_seg)){
        fit1 <- lm(freq_perm[1:t] ~ time[1:t])
        fit2 <- lm(freq_perm[t:n] ~ time[t:n])
        
        rss_perm[t] <- sum(resid(fit1)^2) + sum(resid(fit2)^2)
      }
      
      # compute gain relative to permuted global fit
      fit_all_perm <- lm(freq_perm ~ time)
      rss_all_perm <- sum(resid(fit_all_perm)^2)
      
      gain_perm <- rss_all_perm - rss_perm
      
      perm_gains[i] <- max(gain_perm, na.rm = TRUE)
    }
    
    # p-value = how often permuted >= observed
    p_value <- mean(perm_gains >= best_gain)
  }
  
  return(list(
    breakpoint = bp,
    time = time[bp],
    
    # effect sizes
    gain = best_gain,
    rel_gain = rel_gain,
    
    # curves (for plotting/debugging)
    gain_curve = gain,
    rss_curve = rss,
    
    # baseline
    rss_all = rss_all,
    
    # permutation results
    perm_gains = perm_gains,
    p_value = p_value
  ))
}

#The performance of the test as defined above was assessed in the following manner

#Loop over the variants of interest (in this example: 3)
for(v in c("variant_1","variant_2","variant_3")){
    df <- tr[[v]]

    #Vizualize the trajectories for all available variant-population pairs (with relative gain and permutation Pvalues included underneath each plot)
    #and scan over them
    for(pop in colnames(df)){
        plot(0,0,type="n",xlab="generations",ylab="derived allele frequency",main=paste0(v," ",pop),ylim=range(df[,pop]),xlim=c(1,500),
             sub=paste0("Relative gain = ",round(break.t[[v]][[pop]]$rel_gain,2),", P value = ",break.t[[v]][[pop]]$p_value))
        for(i in 1:nrow(df)) lines(1:nrow(df),df[,pop])
        abline(v=break.t[[v]][[pop]]$breakpoint)
        scan()
    }
}

#Next, the breakpoints were estimated for all trajectories of the variants of interest like this
result_variant_1 = find_break_rss(variant_1$date_mean, variant_1$af, do_perm = T)
result_variant_2 = find_break_rss(variant_2$date_mean, variant_2$af, do_perm = T)
result_variant_3 = find_break_rss(variant_3$date_mean, variant_3$af, do_perm = T)



