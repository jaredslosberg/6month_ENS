library(scattermore)

fit_periodic_loess_mod <- function (cds_object, theta.v.name, assay.name, y.idx, span = 0.3, length.out = 200, plot = FALSE, plot_loess = T,
          fig.title = NULL, point.size = 2.1, point.alpha = 0.6, line.size = 0.8, 
          line.alpha = 0.8, color.by = NULL, color.vars = NULL, color.name = NULL, x_lab = "θ", 
          y_lab = "y", hue.colors = NULL, ...) 
{
  
  theta.v <- pData(cds_object)[,theta.v.name]
  y <- assay(cds_object, assay.name)[y.idx,]
  
  if(!is.null(color.by)){
    color.vars <- pData(cds_object)[,color.by]
  }
  
  stopifnot(`theta.v need to be between 0 - 2pi.` = (min(theta.v) >= 
                                                       0) & (max(theta.v) <= 2 * pi), `The length of theta.v and y should be the same.` = length(theta.v) == 
              length(y))
  x <- c(theta.v - 2 * pi, theta.v, theta.v + 2 * pi)
  ss.total <- sum(scale(y, scale = FALSE)^2)
  y <- rep(y, 3)
  loess.o <- loess(y ~ x, span = span, ...)
  fitted.v <- loess.o$fitted[(length(theta.v) + 1):(length(theta.v) * 
                                                      2)]
  residual.v <- loess.o$residuals[(length(theta.v) + 1):(length(theta.v) * 
                                                           2)]
  rsquared <- 1 - sum(residual.v^2)/ss.total
  pred.x <- seq(0, 2 * pi, length.out = length.out)
  pred.y <- predict(loess.o, newdata = data.frame(x = pred.x))
  pred.df <- data.frame(x = pred.x, y = pred.y)
  if (plot) {
    if (is.null(fig.title)) 
      fig.title <- paste0("(n=", length(theta.v), ")")
    if (!is.null(color.vars)) {
      color.vars <- factor(color.vars)
      stopifnot(`Length of theta.v does not match color.vars` = length(color.vars) == 
                  length(theta.v))
      if (!is.null(hue.colors)) {
        stopifnot(`Number of colors does not match nlevels of color.vars` = nlevels(factor(color.vars)) == 
                    length(hue.colors))
        color_scale <- scale_color_manual(values = hue.colors, 
                                          name = color.name, limits = levels(color.vars))
      }
      else {
        color_scale <- scale_color_discrete(name = color.name, 
                                            limits = levels(color.vars))
      }
      tmp.df <- pData(cds) %>% as.data.frame() %>% cbind(.,data.frame(theta = theta.v, y = y))
      p_aes <- aes_string(x = "theta", y = "y", color = color.by)
    }
    else {
      tmp.df <- pData(cds_object) %>% as.data.frame() %>% cbind(.,data.frame(theta = theta.v, y = y))
      p_aes <- aes_string(x = "theta", y = "y", color = color.by)
      color_scale <- NULL
    }
    fig <- ggplot(data = tmp.df) + geom_scattermore(mapping = p_aes, 
                                                    pointsize = point.size, alpha = point.alpha) +
      labs(x = x_lab, y = y_lab, title = fig.title) +
      scale_x_continuous(limits = c(0, 2 * pi), 
                         breaks = c(0, pi/2, pi, 3 * pi/2, 2 * pi),
                         labels = paste0(c(0, 0.5, 1, 1.5, 2), "π"), name = "θ") + tricycle:::.gg_theme
    
    if(plot_loess == T){
      fig <- fig + geom_path(data = pred.df, 
                  aes_string(x = "x", y = "y"), size = line.size, alpha = line.alpha) 
         
    }
    return(list(fitted = fitted.v, residual = residual.v, 
                pred.df = pred.df, loess.o = loess.o, rsquared = rsquared, 
                fig = fig))
  }
  return(list(fitted = fitted.v, residual = residual.v, pred.df = pred.df, 
              loess.o = loess.o, rsquared = rsquared))
}