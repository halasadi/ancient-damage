

###########   fixed factor FLASH:  on Lindo data  #########################

signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))

signature_set <- colnames(signature_counts)

topic_fit <- maptpx::topics(signature_counts, K=2);

fit_factors <- topic_fit$theta;
fit_omega <- topic_fit$omega;

library(flashr)


res <- signature_counts

K <- 2
theta_mat <- as.numeric()

for(k in 1:K){
    flash_fit <- flash(t(res),
                       tol=1e-5, maxiter_r1 = 100,
                       partype="constant",
                       factor_value = fit_omega[,k], fix_factor = TRUE,
                       nonnegative=FALSE,
                       ash_para = list(control=list(maxiter=1000)))

     lfit <- flash_fit$l;
     theta_fit <- lfit/sum(lfit);
     theta_mat <- cbind(theta_mat, theta_fit)

     res <- res - fit_omega[,k]%*%t(theta_fit)
}

plot(fit_factors[,1], col="red", type="l")
lines(theta_mat[,1], col="blue")

plot(fit_factors[,2], col="red", type="l")
lines(theta_mat[,2], col="blue")


