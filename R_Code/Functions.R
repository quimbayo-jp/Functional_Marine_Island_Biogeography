# Functions ---------------------------------------------------------------------------------------- 
# This code was created all the function used during the analysis 
# and production of the figures

# This function was created to rescale all prectidors considering within of the models --------------
rescale_variables <- function (x) {(x - mean(x))/(1*sd (x))}

# This function was created to examine correlation among the predictors -----------------------------
# Add Histograms in function Pairs
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
# Add Coefficients of correlation
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1)
}

# This function was created to explore the distribution of each variable------------------------------
distribution <- function(x, FIndex){
  par (mfrow=c(1,2))
  hist(x, prob=T, xlab = FIndex, main = NA)
  lines(density(x), col="red")
  
  hist (log (x+1), prob=T, xlab=FIndex, main = NA)
  lines(density(x), col="red")
}

# This function was created to represent the results from the models ---------------------------------
figure_models_current <- function (Model,Title,Values,Axis.Inf){
  plot_model(Model, vline.color = "grey98", show.values =T,
             show.p = T, value.offset = .3, transform = NULL,
             title = Title,
             axis.labels = c("Iso", "Cur.Area", "Age")) +
    #geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=1.5)+
    scale_y_continuous(limits = Values) +
    theme_light() +
    theme (axis.text.x=element_text(size=13, angle=0, family = "sans"),
           axis.text.y = Axis.Inf,
           axis.title.y = element_text(size = rel(1.5), angle = 90, family = "sans"),
           axis.title.x = element_blank())
}

figure_models_past <- function (Model, Values, Axis.Inf){
  plot_model(Model, vline.color = "grey98", show.values =T,
             show.p = T, value.offset = .3, transform = NULL,
             axis.labels = c("Iso", "Past.Area","Age")) +
    #geom_vline(xintercept = 0, linetype="dotted", color = "blue", size=1.5)+
    scale_y_continuous(limits = Values) +
    theme_light() +
    theme (axis.text.x=element_text(size=13, angle=0, family = "sans"),
           axis.text.y = Axis.Inf,
           axis.title.y = element_text(size = rel(1.5), angle = 90, family = "sans"),
           axis.title.x = element_blank())
}



# This function was created to create the maps of functional index ---------------------------------
map_relative_richness <- function (db, column_interest, n_title, name_scale){
  ggplot () +
    geom_polygon(data=ocean_poly_proj, aes(long,lat, group=group), color="black", fill="grey90", size=0.25) +
    geom_polygon(data=land_poly_proj, aes(long,lat, group=group), color="black", fill="white", size=0.25) +
    geom_point(data=as.data.frame(db), aes(fill=column_interest, X.prj, Y.prj),size=3, shape=21)+
    labs(fill=name_scale)+ # \n this code divide the text in two lines
    scale_fill_distiller(palette = "RdYlBu")+
    ggtitle(n_title) +
    theme (axis.text.x=element_blank(),
           axis.text.y = element_blank(),
           axis.title.y = element_blank(),
           axis.title.x = element_blank(),
           axis.ticks = element_blank(),
           panel.background = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.title=element_blank(),
           legend.position="bottom",
           legend.key=element_rect(fill = "white"),
           legend.key.height = unit(0.8,"cm"), # Size of the choropleth
           legend.key.size = unit(1.3, "cm"),
           legend.text=element_text(colour = "black", size = 13, family = "sans"), 
           plot.title =element_text(colour = "black", size = 15, family = "sans", hjust = 0.5))
  
}

# This function was created to plot the linear correlations between species richness and functional indices ------------
Lin_Correlations <- function (db, xvalue, yvalue, categoric.int, method_curve, formula_curve, name.legend, break.scale, 
                      label.legend, color.pch, xlabel, ylabel, axisname1, axisname2, position.text.y, pos.leng) {
  ggplot (db, aes(x=xvalue, y=yvalue))+
    geom_point(size=3, aes(color=categoric.int, shape=categoric.int)) +
    scale_shape_manual(values = c(15,16,17,18))+
    geom_smooth(se = T, method =method_curve,  formula = formula_curve)+
    stat_poly_eq(position = "identity",formula = y ~ splines::bs(x, 3), 
                 aes(label = paste(..rr.label.., stat(p.value.label), sep = "~~~")), 
                 parse = TRUE, label.x = "right", label.y = position.text.y, vstep=0.06, hjust = 1, coef.digits = 2) +
    scale_colour_manual (name=name.legend, breaks = break.scale, labels = label.legend, values=color.pch) +
    labs (x=xlabel, y=ylabel) +
    theme_light() +
    theme (axis.text.y  = element_text(size=10, angle=0, family = "sans"),
           axis.title.y  = axisname1,
           axis.text.x   = element_text(size=10, angle=0, family = "sans"),
           axis.title.x  = axisname2,
           legend.position = pos.leng)
  
}


