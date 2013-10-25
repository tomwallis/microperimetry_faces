# blur slope differences from zero:

source(paste0(getwd(),'/funs/setup_for_calc.R'))

# Are blur slopes different to zero?

# is group 1 (CS) slope different to zero?
(cs.mid <- mean(params$mu[,1,2]))
(cs.ci <- hdi(params$mu[,1,2]))

# is group 2 (LV) slope different to zero?
(lv.mid <- mean(params$mu[,2,2]))
(lv.ci <- hdi(params$mu[,2,2]))


# is group 3 (LV:F) slope different to zero?
(f.mid <- mean(params$mu[,3,2]))
(f.ci <- hdi(params$mu[,3,2]))