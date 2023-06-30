
library(spacrt)
library(magrittr)

n <- 50; p <- 3; B <- 2000; normalize <- FALSE; return_resamples <- FALSE
X_on_Z_fam <- "binomial"
Y_on_Z_fam <- "poisson"

results.dCRT <- list()
results.spaCRT <- list()

system.time({
  for(k in 1:500){
    set.seed(k)

    data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
                 Y = rpois(n = n, lambda = 1),
                 Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))

    results.dCRT[[k]] <- spacrt::dCRT(data, X_on_Z_fam, Y_on_Z_fam, B,
                                      normalize, return_resamples)
    results.spaCRT[[k]] <- spacrt::spaCRT(data, X_on_Z_fam, Y_on_Z_fam, normalize)

    print(k)
  }
})

df.dCRT <- data.frame(x = results.dCRT %>%
                            lapply(function(val) val$p_value) %>% unlist())
df.spaCRT <- data.frame(x = results.spaCRT %>%
                              lapply(function(val) val$p_value) %>% unlist())

(hist.dCRT <- df.dCRT %>%
                ggplot(aes(x=x)) +
                geom_histogram(bins = 10, fill="#69b3a2",
                               color="#e9ecef", alpha=0.7) +
                xlab('p-values') +
                ylab('Count') +
                ylim(c(0,70)) +
                ggtitle("dCRT") +
                theme_light() +
                theme(plot.title = element_text(hjust = 0.5, face = 'bold')))

(hist.spaCRT <- df.spaCRT %>%
                  ggplot(aes(x=x)) +
                  geom_histogram(bins = 10, fill="#69b3a2",
                                 color="#e9ecef", alpha=0.7) +
                  xlab('p-values') +
                  ylab('Count') +
                  ylim(c(0,70)) +
                  ggtitle("spaCRT") +
                  theme_light() +
                  theme(plot.title = element_text(hjust = 0.5, face = 'bold')))

gridExtra::grid.arrange(hist.dCRT, hist.spaCRT, ncol = 2)

