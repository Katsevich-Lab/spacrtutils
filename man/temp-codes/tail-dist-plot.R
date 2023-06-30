

n <- 40; p <- 3; B <- 30000; normalize <- FALSE; return_resamples <- FALSE
X_on_Z_fam <- "binomial"
Y_on_Z_fam <- "poisson"

data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
             Y = rpois(n = n, lambda = 1),
             Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))

X <- data$X; Y <- data$Y; Z <- data$Z
n <- length(X)

# fit X on Z and Y on Z regressions
X_on_Z_fit <- stats::glm(X ~ Z, family = X_on_Z_fam)
Y_on_Z_fit <- stats::glm(Y ~ Z, family = Y_on_Z_fam)

prod_resids <- (X - X_on_Z_fit$fitted.values)*(Y - Y_on_Z_fit$fitted.values)
# compute the test statistic
test_stat <- 1/sqrt(n) * sum(prod_resids)

prod_resid_resamp <- c()

for(b in 1:B){
  # resampling X from X|Z
  resamp_X <- stats::rbinom(n = n, size = 1, prob = X_on_Z_fit$fitted.values)
  # compute the products of residuals for each resampled observation
  prod_resid_resamp[b] <- 1/sqrt(n) * sum((resamp_X - X_on_Z_fit$fitted.values)*
                                            (Y - Y_on_Z_fit$fitted.values))

  print(b)
}

fun.ecdf <- ecdf(prod_resid_resamp)
my.ecdf <- fun.ecdf(sort(prod_resid_resamp))
plot((1 - my.ecdf) ~ sort(prod_resid_resamp), type = 'l')


res.spaCRT <- spaCRT(data, X_on_Z_fam, Y_on_Z_fam, normalize)

sp.approx <- (sort(prod_resid_resamp)) %>% sapply(function(val){
  return(1-res.spaCRT$cdf(val + 1/sqrt(n)*sum(X_on_Z_fit$fitted.values*(Y-Y_on_Z_fit$fitted.values)),
                          P = X_on_Z_fit$fitted.values, W = Y - Y_on_Z_fit$fitted.values))
})

# sp.approx <- (sort(prod_resid_resamp)) %>% sapply(function(val) return(1 - res.spaCRT$cdf(val, X_on_Z_fit$fitted.values, Y - Y_on_Z_fit$fitted.values)))


normal.approx <- 1 - pnorm(sort(prod_resid_resamp),
                           mean = mean(prod_resid_resamp), sd = sd(prod_resid_resamp))



df <- data.frame(prod.res = sort(prod_resid_resamp),
                 actual = (1 - my.ecdf),
                 normal.aprrox = normal.approx,
                 sp.approx = sp.approx)

df[(B*0.9):B, ] %>%
  reshape2::melt(id = 'prod.res') %>%
  ggplot(aes(x = prod.res)) +
  geom_line(aes(y = value, color = variable), linewidth = 0.6) +
  labs(x = "Product Residuals", y = "Upper Tail Prob. (log scale)", color = " ") +
  scale_color_manual(labels = c("Est. True Prob. (Bootstrap)",
                                "Normal Approximation",
                                "SP Approximation"),
                     values = c("red","limegreen","blue")) +
  scale_y_log10() +
  annotation_logticks() +
  theme_light() +
  theme(text = element_text(family="", size=13),
        axis.title.y = element_text(size=15, face="bold",
                                    margin=margin(r=15)),
        axis.title.x = element_text(size=15, face="bold",
                                    margin=margin(t=15)),
        axis.text = element_text(size=11, face="bold", color="black"),
        legend.title = element_text(face="bold"),
        legend.title.align = 0.5,
        legend.key.size = unit(1,"line"),
        legend.position = 'bottom',
        legend.background = element_rect(fill = "white",
                                         linewidth = 1,
                                         linetype = "solid"),
        legend.box.background = element_rect(colour = "black",
                                             linewidth = 1))





















