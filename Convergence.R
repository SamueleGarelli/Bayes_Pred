set.seed(1)

library(DescTools) # needed to compute numerical integrals with the AUC function
library(latex2exp) # needed to write in LaTeX on plots
library(ggplot2) # for plots

# This script calls two functions that have been coded in Rc++ and are defined in the Convergence.cpp file. These are gauss_path_1 and cop_path_1.
# Hence, the file Convergence.cpp has to be sourced in the current R session in order to be executable. 
# There are two different ways to source Rc++ code. Either by running the command Rcpp::sourceCpp("path of the file") or by opening the Rc++ file 
# and clicking on the Source button on the top right-hand corner of the R script (this will just run the aforementioned command). You should use
# the first method if you are working on the R console, since the Source button is available in the R-studio graphical interface only
# (e.g. this might happen if you are working on a remote server). In case you are using R-studio, I recommend the second way. I find it easier
# to avoid mistakes in accessing the file's path.
# Note that, when sourcing the Convergence.cpp files other functions will be loaded into your R environment in addition to the two needed for
# computing paths. Such functions are called inside the predictive resampling algorithm so they are necessary.

# Gaussian predictive

N = 500
grd = seq(from = -4, to = 3, length.out = 1000)
gauss_1 = gauss_path_1(N,grd)

d_p_gauss = c()

for (i in 1:N){
  
  d_p_gauss[i] = AUC(x = grd, y = abs(gauss_1$p[1,] - gauss_1$p[i,]))
  
}

p_1 = ggplot() + 
  geom_line(data = data.frame(x = 1:N,y = d_p_gauss[1:length(d_p_gauss)]), aes(x = x, y = y), linewidth = 0.4) + 
  labs(x = "n", y = TeX("$||f_0 - f_n||_1$"), color = " ") + 
  ggtitle("Gaussian predictive") +
  theme(plot.title = element_text(hjust = 0.5))

p_1

colors_gauss = c("black","red","blue","orange","green","brown","purple")
p_2 = ggplot() +
  geom_line(data = data.frame(x=grd,y=gauss_1$p[1,], group = "n = 1"), aes(x = x, y = y, color = "n = 1"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=gauss_1$p[2,], group = "n = 2"), aes(x = x, y = y, color = "n = 2"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=gauss_1$p[3,], group = "n = 3"), aes(x = x, y = y, color = "n = 3"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=gauss_1$p[4,], group = "n = 4"), aes(x = x, y = y, color = "n = 4"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=gauss_1$p[10,], group = "n = 10"), aes(x = x, y = y, color = "n = 10"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=gauss_1$p[50,], group = "n = 50"), aes(x = x, y = y, color = "n = 50"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=gauss_1$p[500,], group = "n = 500"), aes(x = x, y = y, color = "n = 500"), linewidth = 0.4) +
  ggtitle("Gaussian predictive") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "x", y = TeX("$f_n(x)$"), color = " ") + 
  scale_color_manual(values = colors_gauss,labels = c("n = 1", "n = 2", "n = 3", "n = 4", "n = 10","n = 50", "n = 500")) + 
  theme(legend.position = c(0.8,0.7),
        legend.key.size = unit(0.01, "cm"),
        legend.key.width = unit(0.3,"cm"),
        legend.title = element_blank(), # Remove legend title
        legend.background = element_rect(fill = "white")) # Set legend background color

p_2

# Copula-based predictive

grd = seq(from = -3, to = 3, length.out = 1000)

N = 8000
cop_1 = cop_path_1(N,grd,rho = 0.9) # rho is chosen arbitrarily since we are not doing any predictive resampling or density estimation

d_p_cop = c()

for (i in 1:N){
  
  d_p_cop[i] = AUC(x = grd, y = abs(cop_1$p[1,] - cop_1$p[i,]))
  
}

p_3 = ggplot() + 
  geom_line(data = data.frame(x = 1:N,y = d_p_cop[1:length(d_p_cop)]), aes(x = x, y = y), linewidth = 0.4) + 
  labs(x = "n", y = TeX("$||f_0 - f_n||_1$"), color = " ") + 
  ggtitle(TeX("Copula-based predictive, $\\alpha_n = (2 - 1/n)/(n+1)$")) +
  theme(plot.title = element_text(hjust = 0.5))

p_3

colors_cop = c("black","red","blue","orange","green","brown","purple")
p_4 = ggplot() +
  geom_line(data = data.frame(x=grd,y=cop_1$p[1,], group = "n = 1"), aes(x = x, y = y, color = "n = 1"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=cop_1$p[2,], group = "n = 2"), aes(x = x, y = y, color = "n = 2"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=cop_1$p[4,], group = "n = 4"), aes(x = x, y = y, color = "n = 4"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=cop_1$p[10,], group = "n = 10"), aes(x = x, y = y, color = "n = 10"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=cop_1$p[500,], group = "n = 500"), aes(x = x, y = y, color = "n = 500"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=cop_1$p[3000,], group = "n = 3000"), aes(x = x, y = y, color = "n = 3000"), linewidth = 0.4) +
  geom_line(data = data.frame(x=grd,y=cop_1$p[8000,], group = "n = 8000"), aes(x = x, y = y, color = "n = 8000"), linewidth = 0.4) +
  ggtitle(TeX("Copula-based predictive, $\\alpha_n = (2 - 1/n)/(n+1)$")) +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "x", y = TeX("$f_n(x)$"), color = " ") + 
  scale_color_manual(values = colors_cop,labels = c("n = 1", "n = 2", "n = 4", "n = 10", "n = 500", "n = 3000","n = 8000")) + 
  theme(legend.position = c(0.8,0.7),
        legend.key.size = unit(0.01, "cm"),
        legend.key.width = unit(0.3,"cm"),
        legend.title = element_blank(), # Remove legend title
        legend.background = element_rect(fill = "white")) # Set legend background color

p_4 

# Copula-based predictive (new alpha)
grd = seq(from = -3, to = 3, length.out = 1000)

N = 400
cop_na_1 = cop_path_1_new_alpha(N,grd,rho = 0.9) # rho is chosen arbitrarily since we are not doing any predictive resampling or density estimation

d_p_cop_na = c()

for (i in 1:N){
  
  d_p_cop_na[i] = AUC(x = grd, y = abs(cop_na_1$p[1,] - cop_na_1$p[i,]))
  
}

p_5 = ggplot() + 
        geom_line(data = data.frame(x = 1:N,y = d_p_cop_na[1:length(d_p_cop_na)]), aes(x = x, y = y), linewidth = 0.8) + 
        labs(x = "n", y = TeX("$||f_0 - f_n||_1$"), color = " ") + 
        ggtitle(TeX("Copula-based predictive, $\\alpha_n = (2 - 1/n)/((n+1)(log(n+1))^2)$")) +
        theme(plot.title = element_text(hjust = 0.5))

p_5

colors_cop = c("black","red","blue","orange","green","brown","purple")
p_6 = ggplot() +
        geom_line(data = data.frame(x=grd,y=cop_na_1$p[1,], group = "n = 1"), aes(x = x, y = y, color = "n = 1"), linewidth = 0.8) +
        geom_line(data = data.frame(x=grd,y=cop_na_1$p[5,], group = "n = 5"), aes(x = x, y = y, color = "n = 5"), linewidth = 0.8) +
        geom_line(data = data.frame(x=grd,y=cop_na_1$p[10,], group = "n = 10"), aes(x = x, y = y, color = "n = 10"), linewidth = 0.8) +
        geom_line(data = data.frame(x=grd,y=cop_na_1$p[25,], group = "n = 25"), aes(x = x, y = y, color = "n = 25"), linewidth = 0.8) +
        geom_line(data = data.frame(x=grd,y=cop_na_1$p[50,], group = "n = 50"), aes(x = x, y = y, color = "n = 50"), linewidth = 0.8) +
        geom_line(data = data.frame(x=grd,y=cop_na_1$p[100,], group = "n = 100"), aes(x = x, y = y, color = "n = 100"), linewidth = 0.8) +
        geom_line(data = data.frame(x=grd,y=cop_na_1$p[400,], group = "n = 400"), aes(x = x, y = y, color = "n = 400"), linewidth = 0.8) +
        ggtitle(TeX("Copula-based predictive, $\\alpha_n = (2 - 1/n)/((n+1)(log(n+1))^2)$")) +
        theme(plot.title = element_text(hjust = 0.5)) + 
        labs(x = "x", y = TeX("$f_n(x)$"), color = " ") + 
        scale_color_manual(values = colors_cop,labels = c("n = 1", "n = 5", "n = 10", "n = 25", "n = 50", "n = 100","n = 400")) + 
        theme(legend.position = c(0.8,0.7),    # Move the legend to the bottom
        legend.title = element_blank(), # Remove legend title
        legend.background = element_rect(fill = "white")) # Set legend background color

p_6

